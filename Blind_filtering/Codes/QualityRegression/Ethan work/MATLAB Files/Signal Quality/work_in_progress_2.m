function [classifyResult] = challenge(recordName)
%%
% INPUT:
% recordName: input PCG signal name
%
% OUTPUT:
% classifyResult: Classification result(s)
%
% Contact:
% Morteza Zabihi (morteza.zabihi@gmail.com) && Ali Bahrami Rad(abahramir@yahoo.com)
% Black Swan Team (April 2016)
% This code is released under the MIT License (MIT) (http://opensource.org/licenses/MIT)
%
%% Parameters--------------------------------------------------------------
fs = 2000;
hamming = @(N)(0.54-0.46*cos(2*pi*[0:N-1].'/(N-1)));

load('ANNcinc.mat')
rng(100);
TF = strcmp('.wav',recordName(end-3:end));
if TF ~=1
    recordName = [recordName '.wav'];
end

%% Feature Extraction------------------------------------------------------,
[x, ~] = audioread(recordName);  % load data
% ---------------------------------------------------------------------
[feature] = feature_extraction(x, hamming);% 23 features
end


function [feature] = feature_extraction(x, hamming)
%% 
% INPUT:
% x: Raw data
% hamming: hamming window
%
% OUTPUT:
% feature: 18 features extracted from time-frequnecy domain
%
% Contact:
% Morteza Zabihi (morteza.zabihi@gmail.com) && Ali Bahrami Rad(abahramir@yahoo.com)
% Black Swan Team (July 2016)
% This code is released under the MIT License (MIT) (http://opensource.org/licenses/MIT)
%
%%
Fs = 2000;

%% Linear predictive coefficient 
% Linear Predictive Coefficient (LPC): 10th order linear predictor with all coefficients are used as features
[lp,~] = lpc(double(x),10);

%% Entropy based features
nbins = 50;
[counts,centers] = hist(x,nbins);
thr = mean((diff(centers))/2);
centers(2,:) = centers + thr;
for i=1:length(x)
    p1 = find(x(i)<=centers(2,:));
    if ~isempty(p1)
        p(i) = counts(p1(1));
    else
        p(i) = counts(end);
    end
end
p = p / length(x);
q = 2;
Shannon = -1*sum(p.*log(p));     
Tsallis = (1/(q-1))*sum(1 - p.^q); 

%% MFCC
[MFCCs, ~, ~ ] = mfcc( x, Fs, 0.025*Fs, 0.010*Fs, 0.97, hamming, [1 1000], 26, 14, 22 );

min_MFCC = mean(min(MFCCs)); 
m4 = moment((max(MFCCs)),4,2);       
max_MFCC = abs(m4.^(1/4));  
m4 = moment((skewness(MFCCs)),4,2);  
skew_MFCC = abs(m4.^(1/4));   
                                                             
%% Wavelet transform based feature
level = 5;
[c,l] = wavedec(x,level,'db4');

a5 = appcoef(c,l,'db4');
det5 = detcoef(c,l,level);
det4 = detcoef(c,l,level-1);
det3 = detcoef(c,l,level-2);

wav_log = log2(var(det3)); 
%-------------------------------------------------------------
[counts1,centers1] = hist(a5,nbins);
thr = mean((diff(centers1))/2);
centers1(2,:) = centers1 + thr;
for i=1:length(x)
    p1 = find(x(i)<=centers1(2,:));
    if ~isempty(p1)
        p_1(i) = counts1(p1(1));
    else
        p_1(i) = counts1(end);
    end
end
p_1 = p_1 / length(x);
Shannon_a5 = -1*sum(p_1.*log(p_1));                             
%-------------------------------------------------------------
[counts2,centers2] = hist(det5,nbins);
thr2 = mean((diff(centers2))/2);
centers2(2,:) = centers2 + thr2;
for i=1:length(x)
    p1 = find(x(i)<=centers2(2,:));
    if ~isempty(p1)
        p2(i) = counts2(p1(1));
    else
        p2(i) = counts2(end);
    end
end
p2 = p2 / length(x);
q = 2;
Renyi_d5 = (1/(q-1))*log(sum(p2.^q));                        
%-------------------------------------------------------------´
[counts3,centers3] = hist(det4,nbins);
thr3 = mean((diff(centers3))/2);
centers3(2,:) = centers3 + thr3;
for i=1:length(x)
    p1 = find(x(i)<=centers3(2,:));
    if ~isempty(p1)
        p3(i) = counts3(p1(1));
    else
        p3(i) = counts3(end);
    end
end
p3 = p3 / length(x);
Shannon_d4 = -1*sum(p3.*log(p3));                             
%% ------------------------------------------------------------------------
%% Features extracted over power spectral density
xx = downsample(x, 16);
% the areas under the curve of the frequency domain
[Pfb1, Pfb2, Pfb3, Pfb4, Pfb5, Pfb6,Pfb7, Pfb8, Pfb9, Pfb10] = Spectral(xx);                               
                     
N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(2*pi*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);

f = Fs*(0:(N/2))/N;
% modified power spectral densit centroid.
frequency_centroid = (sum(f'.*(psdx.^2)))/(sum(psdx.^2));   



feature = [f1 f5 f6 f7 f14 f15 f161821232425 frequency_centroid,...
            Shannon Tsallis Shannon_a5 Renyi_d5 Shannon_d4];
        
end

%
function [Pfb1, Pfb2, Pfb3, Pfb4, Pfb5, Pfb6,Pfb7, Pfb8, Pfb9, Pfb10] = Spectral(x)
%%
% INPUT:
% x:Raw data
%
% OUTPUT:
% Pfb8, Pfb9, Pfb10: 3 features extracted from frequency domain
%
% Contact:
% Morteza Zabihi (morteza.zabihi@gmail.com) && Ali Bahrami Rad(abahramir@yahoo.com)
% Black Swan Team (April 2016)
% This code is released under the MIT License (MIT) (http://opensource.org/licenses/MIT)
%
%%
fb1 =  [0 0.1];
fb2 =  [0.1 0.2];
fb3 =  [0.2 0.3];
fb4 =  [0.3 0.4];
fb5 =  [0.4 0.5];
fb6 =  [0.5 0.6];
fb7 =  [0.6 0.7];
fb8 =  [0.7 0.8];
fb9 =  [0.8 0.9];
fb10 = [0.9 1];
%%
% uses a hamming window
[PSD,F] = pwelch(x); 
% find the indexes corresponding bands 
ifb1 = (F>=fb1(1)) & (F<=fb1(2));
ifb2 = (F>=fb2(1)) & (F<=fb2(2));
ifb3 = (F>=fb3(1)) & (F<=fb3(2));
ifb4 = (F>=fb4(1)) & (F<=fb4(2));
ifb5 = (F>=fb5(1)) & (F<=fb5(2));
ifb6 = (F>=fb6(1)) & (F<=fb6(2)); 
ifb7 = (F>=fb7(1)) & (F<=fb7(2));

ifb8 =  (F>=fb8(1))  & (F<=fb8(2));
ifb9 =  (F>=fb9(1))  & (F<=fb9(2));
ifb10 = (F>=fb10(1)) & (F<=fb10(2));

%calculate areas, within the freq bands (ms^2) 
Ifb1 = trapz(F(ifb1),PSD(ifb1));
Ifb2 = trapz(F(ifb2),PSD(ifb2));
Ifb3 = trapz(F(ifb3),PSD(ifb3));
Ifb4 = trapz(F(ifb4),PSD(ifb4));
Ifb5 = trapz(F(ifb5),PSD(ifb5));
Ifb6 = trapz(F(ifb6),PSD(ifb6));
Ifb7 = trapz(F(ifb7),PSD(ifb7));

Ifb8  = trapz(F(ifb8),PSD(ifb8));
Ifb9  = trapz(F(ifb9),PSD(ifb9));
Ifb10 = trapz(F(ifb10),PSD(ifb10));

aTotal = Ifb1+Ifb2+Ifb3+Ifb4+Ifb5+Ifb6+Ifb7+Ifb8+Ifb9+Ifb10;
% calculate areas relative to the total area (%)
Pfb1 =(Ifb1/aTotal)*100;
Pfb2 =(Ifb2/aTotal)*100;
Pfb3 =(Ifb3/aTotal)*100; 
Pfb4 =(Ifb4/aTotal)*100;
Pfb5 =(Ifb5/aTotal)*100;
Pfb6 =(Ifb6/aTotal)*100;
Pfb7 =(Ifb7/aTotal)*100;
Pfb8 =(Ifb8/aTotal)*100;
Pfb9 =(Ifb9/aTotal)*100;
Pfb10 =(Ifb10/aTotal)*100;
end