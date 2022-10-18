%%%%%%%%%%%%%%%%%%%
%% Initializing
clear
close all
clc
DA='..\10sec\' %%% directory address
files=dir('*.mp3'); %% Defining address of the folder containing the audio recorde



fcl=100;  %% Low cutoff frequency
fch=1000;    %% High cutoff frequency
OutputSheet=cell(size(files,1)+1,7);
OutputSheet(1,1:16)={'name', 'mfcc1 mean', 'mfcc2 mean', 'mfcc3 mean', 'mfcc4 mean', 'mfcc5 mean', 'mfcc6 mean', 'mfcc7 mean', 'mfcc8 mean', 'mfcc9 mean', 'mfcc10 mean', 'mfcc11 mean', 'mfcc12 mean', 'mfcc13 mean', 'mfcc14 mean', 'mfcc15 mean'};
%% Doing for each of the records ... 
for i=1:size(files,1)
    [data,fs]=audioread([DA files(i).name]);    %% read the record
    
    %% MFCC 
    mfcc_mean = mean(mfcc(data,fs, 'NumCoeffs', 15,'LogEnergy','Ignore'));
   
    %% Feature extracting
    OutputSheet(i+1,1)=cellstr(files(i).name);
    OutputSheet(i+1,2:16)=num2cell(mfcc_mean);

end

writecell(OutputSheet,'..\results\MFCC15 features.xls');
