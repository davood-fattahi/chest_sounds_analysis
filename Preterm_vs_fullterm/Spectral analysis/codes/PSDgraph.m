%%%%%%%%%%%%%%%%%%%
% this function provide and compare the psd graphs of preterm vs fullterm
% cases.
% 
% Davood Fattahi, Feb 2020

%% Initializing
clear
close all
clc

% frequency range 
fcl=100;
fch=2000;
% DA='..\Samples\';
DA='..\10sec\';


files=dir([DA '*.mp3']); %% Defining address of the folder containing the audio records

%%%% Pre-allocation ...
%P=zeros(1,size(files,1));
name=cell(size(files,1),1);

for i=1:size(files,1)
    %%% Display the progress percentage 
    clc;disp([num2str(i) '/' num2str(size(files,1)) ', ' num2str(floor(100*i/size(files,1))) ' %']) 
    [data,fs]=audioread([DA files(i).name]);    %% read the record
%     [data,fs]=audioread(['10sec\' files(i).name]);    %% read the record
    name{i}=files(i).name;
    
    %%% find non-clip intervals
%   data=nonclip(data,fs,1,0.1,.20,'concat');
    
%      %%
%     % Bandpass filtering 1 :     
%     dlp=fdesign.lowpass('Fp,Fst,Ap,Ast', .9*fch, 1.1*fch, .1, 60,fs);
%     Hdlp = design(dlp,'equiripple');
%     % fvtool(Hdlp)
%     % measure(Hd)
% 
%     dhp = fdesign.highpass('Fst,Fp,Ast,Ap', .9*fcl, 1.1*fcl, 60, .1, fs);
%     Hdhp = design(dhp,'equiripple');
%     % fvtool(Hdhp)
%     % measure(Hdhp)
% 
%     data=filter(Hdlp,(filter(Hdhp,data)));
%     
%     
%     %% Normalizing the data
%     data=data./((std(data,0,1)));
%     
    
    %% pwelch
    %%% parameter setting for Welch
    wind = ones(1,floor(0.5*fs)); % 500 m sec
    nover = floor(length(wind)*.5); % 0.5 overlap
    nfft = 2^(nextpow2(length(wind))-1)/2; % nfft

    %%% Power Spectral Density
   
    [pxx,f]=pwelch(data,wind,nover,nfft,fs); % Welch
    P(:,i) = pxx./bandpower(pxx,f,[fcl fch],'psd');
end
% frequency range for plotting
ff=f(f<=fch);
P=P(f<=fch,:);

[Lnum1,Ltxt1,Lraw1]=xlsread('..\Demographics1.xlsx','Sheet1'); %% Pre terms
[Lnum2,Ltxt2,Lraw2]=xlsread('..\Demographics1.xlsx','Sheet2');  %% full terms

Lname1=Ltxt1(2:end,1);  %% Pre terms
Lname2=Ltxt2(2:end,1);  %% full terms
Lnum=[Lnum1;Lnum2];   %% all terms
Lname=[Lname1;Lname2];  %% all terms

%% removing extensions form file names
for i=1:size(name,1)
    if name{i}(end-3:end-2)=='.m' 
       name{i}(end-3:end)=[];
    end
end


%% Reordering the features according the labels

%finding common cases in labels and features - pre
[ind1,loc1]=ismember(Lname1,name);

% Removing the labels not existing in the features - pre
L1=Lnum1(loc1~=0,:);

% Removing the features not existing in the labels - pre
P1=P(:,loc1(loc1~=0));
name1=name(loc1(loc1~=0),:);

%finding common cases in labels and features - full
[ind2,loc2]=ismember(Lname2,name);

% Removing the labels not existing in the features - full
L2=Lnum2(loc2~=0,:);

% Removing the features not existing in the labels - full
P2=P(:,loc2(loc2~=0));
name2=name(loc2(loc2~=0),:);

%% Plot the results
%%% power median and confidence interval:
Pav1=median(P1,2);
Pav2=median(P2,2);
dev1=[prctile(P1',75)'-Pav1 Pav1-prctile(P1',25)'];
dev2=[prctile(P2',75)'-Pav2 Pav2-prctile(P2',25)'];

figure; 
shadedErrorBar(ff,Pav1',dev1','lineprops','-b','transparent',1);
hold on
shadedErrorBar(ff,Pav2',dev2','lineprops','-r','transparent',1);
hold off
title('Power Spectral Density');
legend('pre term', 'full term')
xlabel('frequency (Hz)'); ylabel('Normalized PSD')
xlim([100 2000])
xticks(100:200:2000)
saveas(gcf,'..\results\med_psd_intqrtl','epsc')

%%% log of median:
figure; 
plot(ff,10*log10(Pav1)); hold on
plot(ff,10*log10(Pav2));
title('Logarithm of Median Power Spectral Density');
legend('pre term', 'full term')
xlabel('frequency (Hz)'); ylabel('Normalized PSD (dB/Hz)')
xlim([100 2000])
xticks(100:200:2000)
saveas(gcf,'..\results\log_med_psd','epsc')

%%% median of logs:
P1=10*log10(P1);
P2=10*log10(P2);

Pav1=median(P1,2);
Pav2=median(P2,2);
dev1=[prctile(P1',75)'-Pav1 Pav1-prctile(P1',25)'];
dev2=[prctile(P2',75)'-Pav2 Pav2-prctile(P2',25)'];


figure; 
shadedErrorBar(ff,Pav1',dev1','lineprops','-b','transparent',1);
hold on
shadedErrorBar(ff,Pav2',dev2','lineprops','-r','transparent',1);
hold off
title('Power Spectral Density');
legend('pre term', 'full term')
xlabel('frequency (Hz)'); ylabel('Normalized PSD (dB/Hz)')
xlim([100 2000])
xticks(100:200:2000)
saveas(gcf,'..\results\log_med_psd_intqrtl','epsc')





