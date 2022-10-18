%%%%%%%%%%%%%%%%%%%
%% Initializing
clear
% close all
clc
delete 'Results\FreqFixHeart\*.wav'
delete 'Results\FreqFixLung\*.wav'

fs=4000;
Hfr=[50 250];
Lfr=[100 1000];

%% Filter design (do if all the sampling freq is equal, else do it in the loop)
%%% Bandpass filtering for Heart :

[B_H,A_H] = butter(5,[2*Hfr(1)/fs 2*Hfr(2)/fs],'bandpass');
fvtool(B_H, A_H)
hold on
[B_L,A_L] = butter(5,[2*Lfr(1)/fs 2*Lfr(2)/fs],'bandpass');
fvtool(B_L, A_L)

%% Defining address of the folder containing the audio records
DA='..\..\..\chestRecords\currentDataset\';
dataSheet=readcell([DA '00dataSheet.xlsx']);
filesName=dataSheet([false; cell2mat(dataSheet(2:end,5))==1],1);
%% do for each of the records
for i=1:size(filesName,1)
    %%% Display the progress percentage 
    clc;disp([num2str(i) '/' num2str(size(filesName,1)) ', ' num2str(floor(100*i/size(filesName,1))) ' %']) 
    
    %%% Loading the signal
    [data,fs]=audioread([DA filesName{i}]);    %% read the record
  
    %%% Preprocessing - Downsampling
    data=downsample(data,fs./4000); fs=4000;
    data=data(1:10*fs);
    
    %%% Correcting the name - removing extension
    name{i}=filesName{i};
    if strcmp(name{i}(end-3:end-2),'.m') || strcmp(name{i}(end-3:end-2),'.w')
        name{i}(end-3:end)=[];
    end

  
    %% Bandpass filtering for Heart:

    H=filtfilt(B_H,A_H,data);
    H=H./(6.*std(H)); %% normalising by k*variance, to avoid clipping or low amlitude
    HH(i,:)=H;
    %% Bandpass filtering for Lung :

    L=filtfilt(B_L,A_L,data);
    L=L./(15.*std(L)); %% normalising by k*variance, to avoid clipping or low amlitude
    LL(i,:)=L;
% % plot 
%     figure
%     subplot(311)
%     plot(data)
%     subplot(312)
%     plot(H)
%     subplot(313)
%     plot(L)
%% save
    filename = ['Results\FreqFixHeart\' name{i} '_FreqFilt_Heart' '.wav'];
    audiowrite(filename,H,fs)
    filename = ['Results\FreqFixLung\' name{i} '_FreqFilt_Lung' '.wav'];
    audiowrite(filename,L,fs)
end
save('Results\FreqFiltData.mat','HH','LL')




