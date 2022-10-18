%% Zero-detection algorithm
[meanCoefS1, meanCoefSystole, meanCoefS2, meanCoefDiastole, PCG_Hi_Hi] = falkowyFeature(PCG_resampled,assigned_states);
PCG_resampled=PCG_Hi_Hi;
[ qrs_pos,  PCG_Hi_Hi_seg2]=checkingZero2(PCG_Hi_Hi);
[transitCounterNorm,RMSSDall,reducedTransitCounter]=detectNoisySignals2(PCG_Hi_Hi_seg2, qrs_pos);


%% Classification
[ features2 ] = getNewFeatures(assigned_states,PCG_resampled);
decision = newClassifyResults_ALL(features2);


function [PCG_Hi_Hi] = falkowyFeature(PCG)
% input data:
% A) PCG signal

% data processing and output:
% A) calculation of the wavelet transform with the nucleus of Daubechies-2

PCG = butterworth_high_pass_filter(PCG,2,20,1000);
% the beginning of the wavelet transform
falka = dbwavf('db2');
falkaODW =-( (-1).^(1:length(falka)) ).*falka;
PCG_Lo = conv(PCG,falka);
% cut off the tip
PCG_Lo = PCG_Lo(1:end-length(falka)+1);
% downlample
PCG_Lo = PCG_Lo(1:2:end);
PCG_Hi_Hi = conv(PCG_Lo(1:end-length(falka)+1),falkaODW);
% we get rid of the basic words from the plexus
PCG_Hi_Hi = PCG_Hi_Hi(1:end-length(falka)+1);

PCG_Hi_Hi = interp(PCG_Hi_Hi,2);
end



function [qrs_pos,PCG_resampledNew]=checkingZero2(PCG)

len=length(PCG);
PCG_resampledNew=[];

%normalisation
windowN=25;
for j=1:windowN
    PCG_resampledtemp=PCG(1+floor(len/windowN)*(j-1):j*floor(len/windowN));
    if( (max(PCG_resampledtemp)-min(PCG_resampledtemp)) ~=0)
        PCG_resampledNew=[PCG_resampledNew ;...
            (PCG_resampledtemp-mean(PCG_resampledtemp))/(max(PCG_resampledtemp)-min(PCG_resampledtemp))];
    end
end
[mdfint,PCG_resampledNew] =qrs_detect2_PCG(PCG_resampledNew);
peaks=GetPeaks3(mdfint);
qrs_pos=peaks;
end




function  [mdfint2, PCG] =qrs_detect2_PCG(PCG)
%%% I assume that resampling to 1000 is already done
fs = 1000;

[a,b] = size(PCG);
if(a>b)
elseif(b>a)
    PCG=PCG';
end

% == constants
MED_SMOOTH_NB_COEFF = round(fs/100);
% length is 7 for fs=256Hz
INT_NB_COEFF = round(28*fs/256);


PCG = butterworth_low_pass_filter(PCG,2,80,1000, false);
PCG = butterworth_high_pass_filter(PCG,2,30,1000);
bpfPCG=PCG';
bpfPCGNew=[];
len=length(bpfPCG);
windowN=40;
for j=1:windowN
    bpfPCGNewtemp=bpfPCG(1+floor(len/windowN)*(j-1):j*floor(len/windowN));
    bpfPCGNew=[bpfPCGNew  (bpfPCGNewtemp-mean(bpfPCGNewtemp))/(max(bpfPCGNewtemp)-min(bpfPCGNewtemp))];
end

bpfPCG=bpfPCGNew;
% (5) square PCG
sqrPCG = bpfPCG'.*bpfPCG';
% (6) integrate
intPCG = filter(ones(1,INT_NB_COEFF),1,sqrPCG);
% (7) smooth
mdfint = medfilt1(intPCG,MED_SMOOTH_NB_COEFF);
delay  = ceil(INT_NB_COEFF/2);
% remove filter delay for scanning back through PCG
mdfint = circshift(mdfint,-delay);

mdfint2=movingAvgForwBack(mdfint, 40);
end






function peaks=GetPeaks3(sigIR)
flag=0;
tabpeak=zeros(3,1);
slopecount=0;
lastBeatTime=0;
peaks=[];

for i=1:length(sigIR)
    tabpeak(1)=tabpeak(2);
    tabpeak(2)=tabpeak(3);
    tabpeak(3)=sigIR(i);
    if(tabpeak(2)>tabpeak(1) && tabpeak(2)>tabpeak(3) && flag==0)
        flag=1;
        % add window_num * window length -1
        startB=i; 
        startV=tabpeak(2);
    end
    if(flag==1 && tabpeak(3)<tabpeak(2))
        slopecount=slopecount+1;
    end
    
    if(flag==1 && tabpeak(3)>=tabpeak(2))
        if( slopecount>=78 && abs(startV-sigIR(i))>1.5 && startB-lastBeatTime>200)
            peaks=[peaks ; startB];
            lastBeatTime=startB; % +inxMaxSigLocal;
        end
        flag=0;
        slopecount=0;
    end
end

end





function [transitCounterNorm,RMSSDall,reducedTransitCounter]=detectNoisySignals2(PCG_resampled, qrs_pos)

% Number of times which signal intersected the horizontal line
% Horizontal line determined by 0.85 quantile of values of that signal
% divided by teh signal length
% Horizontal line value
quantSigValue=quantile(PCG_resampled,0.85);
% Number of times intersected
transitPosition=find((PCG_resampled(1:end-1)>quantSigValue && PCG_resampled(2:end)<=quantSigValue) ...
    || (PCG_resampled(1:end-1)<quantSigValue && PCG_resampled(2:end)>=quantSigValue));
transitCounter=length(transitPosition);
% Divided by signal length
% Lower than 0.06 is good quality
transitCounterNorm=transitCounter/length(PCG_resampled);

% Root mean square of successive differences
% Lower than 0.026 is good quality
RMSSDall=sqrt(sum(diff(PCG_resampled).^2)/(length(PCG_resampled)-1));

% Value of the 0.58 percentile of signal and the number of detected peaks
% was extracted from the whole number of intersections
%find maximum peak value near the peack detected point
qrs_val=zeros(length(qrs_pos),1);
inxToReject=[];

for i=1:length(qrs_val)
    if(qrs_pos(i)+120<length(PCG_resampled) &&  qrs_pos(i)-120>0)
        [peakValue,peakValueInx]=max(PCG_resampled(qrs_pos(i)-120:qrs_pos(i)+120));
        qrs_val(i)=peakValue;
        inxToReject=[inxToReject  qrs_pos(i)-120+peakValueInx-100:qrs_pos(i)-120+peakValueInx+100];
    end
end
meanPeakValue=quantile(qrs_val,0.58);
inxToReject(inxToReject>length(PCG_resampled) | inxToReject<=0)=[];
PCG_resampled(inxToReject)=[];
transitPosition=find((PCG_resampled(1:end-1)>meanPeakValue && PCG_resampled(2:end)<=meanPeakValue) ...
    || (PCG_resampled(1:end-1)<meanPeakValue && PCG_resampled(2:end)>=meanPeakValue));

transitCounter=length( transitPosition );
% Lower than 18 is good quality
reducedTransitCounter=transitCounter-length( qrs_pos );


end


function [per_peak_window_original,per_peak_window_modified,per_peak_window_modified2]= determineZero(PCG_resampled,qrs_pos)
%2200ms length moving windows with 25% overlap and checked how many peaks are detected in each window.
%Assigned 1 to each window containing 2-4 peaks and 0 in other cases.
window=2200;
overlap=0.25;
peacksDensity=[];
thdown=2;
thup=4;

% we create the ECG index vector
offset=(1:round(overlap*window):size(PCG_resampled));
% we create a vector [0 1 2 ... window-1] - one window long
id1=(0:window-1)';
% The bsxfun function creates an index table - such that in the n-th column there are subsequent samples from the n-th signal window
idx=bsxfun(@plus, id1, offset);
IDX=[];
for j=1:size(idx,2)
    if(idx(end,j)<length(PCG_resampled))
        IDX=[IDX  idx(:,j)];
    end
end

for i=1:size(IDX,2)
    peacksDensity=[peacksDensity ; sum(ismember(qrs_pos,IDX(:,i)))];
end

peacksDensityBin=peacksDensity>=thdown && peacksDensity<=thup;
% Criterion of 65% of windows with score equal to 1
per_peak_window_original=sum(peacksDensityBin)/length(peacksDensityBin);

%5 and 95 percentile heart rate Term: 120-185, 17+: 60-115
% 2 to 4 for 2200ms window is equivalent to 54.5455 to 109.0909
% For newbirns should be 4-7  109.0909 to 190.9091
peacksDensityBin=peacksDensity>=4 && peacksDensity<=7;
per_peak_window_modified=sum(peacksDensityBin)/length(peacksDensityBin);

%Also consider full possible range of heart rate which is max_HR=220; min_HR=70;
% For newbirns should be 4-7  54.5455 to 218.1818
peacksDensityBin=peacksDensity>=2 && peacksDensity<=8;
per_peak_window_modified2=sum(peacksDensityBin)/length(peacksDensityBin);
end

function [coeffRMSSD, coeffZERO, coeffSD1] = getNewFeatures(assigned_states,PCG)
% find the locations with changed states
indx = find(abs(diff(assigned_states))>0);

% for some recordings, there are state zeros at the beginning of assigned_states
if assigned_states(1)>0
    switch assigned_states(1)
        case 4
            K=1;
        case 3
            K=2;
        case 2
            K=3;
        case 1
            K=4;
    end
else
    switch assigned_states(indx(1)+1)
        case 4
            K=1;
        case 3
            K=2;
        case 2
            K=3;
        case 1
            K=0;
    end
    K=K+1;
end

indx2 = indx(K:end);
rem = mod(length(indx2),4);
indx2(end-rem+1:end) = [];
% A is N*4 matrix, the 4 columns save the beginnings of S1, systole, S2 and diastole in the same heart cycle respectively
A = reshape(indx2,4,length(indx2)/4)';
%% RMSSD divided into the number of samples

RMSSD_sys_norm=zeros(size(A,1)-1,1);
RMSSD_dias_norm=zeros(size(A,1)-1,1);
RMSSD_all_norm=zeros(size(A,1)-1,1);

zero_sys_norm=zeros(size(A,1)-1,1);
zero_dias_norm=zeros(size(A,1)-1,1);
zero_all_norm=zeros(size(A,1)-1,1);

SD1_sys=zeros(size(A,1)-1,1);
SD1_dias=zeros(size(A,1)-1,1);
SD1_all=zeros(size(A,1)-1,1);

for i=1:size(A,1)-1
    %normalization
    PCG(A(i,1):A(i+1,1))=(PCG(A(i,1):A(i+1,1))-mean(PCG(A(i,1):A(i+1,1))))/...
        (max(PCG(A(i,1):A(i+1,1)))-min(PCG(A(i,1):A(i+1,1))));
    
    PCGtemp_sys=PCG(A(i,2):A(i,3));
    nsys=length(PCGtemp_sys);
    PCGtemp_sys_diff=diff(PCGtemp_sys);
    xp=PCGtemp_sys_diff;
    xp(end)=[];
    xm=PCGtemp_sys_diff;
    xm(1)=[];
    % SD1
    SD1_sys(i) = std(xp-xm)/sqrt(2);
    
    PCGtemp_dias=PCG(A(i,4):A(i+1,1));
    ndias=length(PCGtemp_dias);
    PCGtemp_dias_diff=diff(PCGtemp_dias);
    xp=PCGtemp_dias_diff;
    xp(end)=[];
    xm=PCGtemp_dias_diff;
    xm(1)=[];
    %SD1
    SD1_dias(i) = std(xp-xm)/sqrt(2);
    PCGtemp_all=[PCGtemp_sys ; PCGtemp_dias];
    nall=length(PCGtemp_all);
    
    RMSSD_sys_norm(i)=  sqrt(sum(diff(PCGtemp_sys).^2)/(nsys-1));
    RMSSD_dias_norm(i) = sqrt(sum(diff(PCGtemp_dias).^2)/(ndias-1));
    RMSSD_all_norm(i) = sqrt(sum(diff(PCGtemp_all).^2)/(nall-1));
    
    zero_sys_norm(i)=sum( (PCGtemp_sys(1:end-1)>0 && PCGtemp_sys(2:end)<=0)...
        || (PCGtemp_sys(1:end-1)<=0 && PCGtemp_sys(2:end)>0))/nsys;
    zero_dias_norm(i)=sum( (PCGtemp_dias(1:end-1)>0 && PCGtemp_dias(2:end)<=0)...
        || (PCGtemp_dias(1:end-1)<=0 && PCGtemp_dias(2:end)>0))/ndias;
    zero_all_norm(i)=sum( (PCGtemp_all(1:end-1)>0 && PCGtemp_all(2:end)<=0)...
        || (PCGtemp_all(1:end-1)<=0 && PCGtemp_all(2:end)>0))/nall;
    
    
    PCGtemp_all_diff=diff(PCGtemp_all);
    xp=PCGtemp_all_diff;
    xp(end)=[];
    xm=PCGtemp_all_diff;
    xm(1)=[];
    %SD1
    SD1_all(i) = std(xp-xm)/sqrt(2);
end

coeffRMSSD=(RMSSD_sys_norm-RMSSD_dias_norm)./RMSSD_all_norm;
coeffZERO=(zero_sys_norm-zero_dias_norm)./zero_all_norm;

coeffSD1=(SD1_sys-SD1_dias)./SD1_all;
end



function [per_abnormal_rmssd,per_abnormal_zero,per_abnormal_sd1]=newClassifyResults_ALL(coeffRMSSD, coeffZERO, coeffSD1)
%thresholds for 3 selected features
th1=0.8; %rmssd
th2=0.8; %sd1
th4=0.6; %zero

% The percentage of bins where the threshold for
% of the trait is abnormal
%p1=0.8;  %rmssd
%p2=0.7; %sd1
%p4=0.8; %zero
per_abnormal_rmssd=sum(abs(coeffRMSSD)>th1)/length(coeffRMSS);
per_abnormal_sd1=sum(abs(coeffSD1)>th2)/legnth(coeffZERO);
per_abnormal_zero=sum(abs(coeffZERO)>th4)/length(coeffZERO);
end