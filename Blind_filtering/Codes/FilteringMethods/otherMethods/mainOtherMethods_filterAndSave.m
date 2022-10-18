%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%% Initializing
clear
% close all
clc

% 
% delete 'TempResults\*\*.wav'
% delete 'TempResults\*\*.wav'


%% Defining address of the folder containing the audio records
DA='..\..\..\chestRecords\currentDataset\';
dataSheet=readcell([DA '00dataSheet.xlsx']);
filesName=dataSheet([false; cell2mat(dataSheet(2:end,5))==1],1);

%% do for each of the records
for i= 1 :size(filesName,1)

    %%% Display the progress percentage 
    clc;disp([num2str(i) '/' num2str(size(filesName,1)) ', ' num2str(floor(100*i/size(filesName,1))) ' %']) 
    
    %%% Loading the signal
    [data,fs]=audioread([DA  '\' filesName{i}]);    %% read the record
    
    %%% Preprocessing - Downsampling to 4000 Hz
    data=downsample(data,fs./4000); fs=4000;
   
    %%% Preprocessing - normalizing
    data=normalize(data);
    
    %%% setting length to 10 sec
    data=data(1:10*fs);
    
        
    [heartSound,lungSound]=emd_separation(data,4000,'emd','log energy');
    heartSound=heartSound./(6.*(std(heartSound,0,1))); % normalization
    lungSound=lungSound./(15.*(std(lungSound,0,1))); % normalization
    audiowrite(['TempResults\emdHeart\' filesName{i}(1:end-4) '_emd_Heart' '.wav'],heartSound,fs)
    audiowrite(['TempResults\emdLung\' filesName{i}(1:end-4) '_emd_Lung' '.wav'],lungSound,fs)

    [heartSound,lungSound]=emd_separation(data,4000,'eemd','log energy');
    heartSound=heartSound./(6.*(std(heartSound,0,1))); % normalization
    lungSound=lungSound./(15.*(std(lungSound,0,1))); % normalization
    audiowrite(['TempResults\eemdHeart\' filesName{i}(1:end-4) '_eemd_Heart' '.wav'],heartSound,fs)
    audiowrite(['TempResults\eemdLung\' filesName{i}(1:end-4) '_eemd_Lung' '.wav'],lungSound,fs)
   
   [heartSound,lungSound]=emd_separation(data,4000,'ceemd','log energy');
    heartSound=heartSound./(6.*(std(heartSound,0,1))); % normalization
    lungSound=lungSound./(15.*(std(lungSound,0,1))); % normalization
    audiowrite(['TempResults\ceemdHeart\' filesName{i}(1:end-4) '_ceemd_Heart' '.wav'],heartSound,fs)
    audiowrite(['TempResults\ceemdLung\' filesName{i}(1:end-4) '_ceemd_Lung' '.wav'],lungSound,fs)
    

    options_nmf.cf = 'kl';
    options_nmf.sparsity = 0;
    xhat= nmf_cluster_modified(data,fs,[], 'Filtering', 'STFT', [], options_nmf);         
    heartSound=xhat(:,1)./(6.*(std(xhat(:,1),0,1))); % normalization
    lungSound=xhat(:,2)./(15.*(std(xhat(:,2),0,1))); % normalization
    audiowrite(['TempResults\nmfcklHeart\' filesName{i}(1:end-4) '_nmfckl_Heart' '.wav'],heartSound,fs)
    audiowrite(['TempResults\nmfcklLung\' filesName{i}(1:end-4) '_nmfckl_Lung' '.wav'],lungSound,fs)

    
    options_nmf.cf = 'ed';
    options_nmf.sparsity = 0;
    xhat= nmf_cluster_modified(data,fs,[], 'Synthesis', 'STFT', [], options_nmf);         
    heartSound=xhat(:,1)./(6.*(std(xhat(:,1),0,1))); % normalization
    lungSound=xhat(:,2)./(15.*(std(xhat(:,2),0,1))); % normalization
    audiowrite(['TempResults\nmfcGenFiltHeart\' filesName{i}(1:end-4) '_nmfcGenFilt_Heart' '.wav'],heartSound,fs)
    audiowrite(['TempResults\nmfcGenFiltLung\' filesName{i}(1:end-4) '_nmfcGenFilt_Lung' '.wav'],lungSound,fs)
      
    xhat= nmf_cluster_modified(data,fs,[], 'No Mask', 'STFT', [], options_nmf);         
    heartSound=xhat(:,1)./(6.*(std(xhat(:,1),0,1))); % normalization
    lungSound=xhat(:,2)./(15.*(std(xhat(:,2),0,1))); % normalization
    audiowrite(['TempResults\nmfcGenNomaskHeart\' filesName{i}(1:end-4) '_nmfcGenNomask_Heart' '.wav'],heartSound,fs)
    audiowrite(['TempResults\nmfcGenNomaskLung\' filesName{i}(1:end-4) '_nmfcGenNomask_Lung' '.wav'],lungSound,fs)
      
    method.options_nmf.cf = 'kl';
    method.options_nmf.sparsity = 2;
    xhat= nmf_cluster_modified(data,fs,[], 'Filtering', 'STFT', [], options_nmf);         
    heartSound=xhat(:,1)./(6.*(std(xhat(:,1),0,1))); % normalization
    lungSound=xhat(:,2)./(15.*(std(xhat(:,2),0,1))); % normalization
    audiowrite(['TempResults\nmfcSparseHeart\' filesName{i}(1:end-4) '_nmfcSparse_Heart' '.wav'],heartSound,fs)
    audiowrite(['TempResults\nmfcSparseLung\' filesName{i}(1:end-4) '_nmfcSparse_Lung' '.wav'],lungSound,fs)
      
    
    [heartSound,lungSound]=singular_spectrum_analysis(data,4000);
    heartSound=heartSound./(6.*(std(heartSound,0,1))); % normalization
    lungSound=lungSound./(15.*(std(lungSound,0,1))); % normalization
    audiowrite(['TempResults\ssaHeart\' filesName{i}(1:end-4) '_ssa_Heart' '.wav'],heartSound,fs)
    audiowrite(['TempResults\ssaLung\' filesName{i}(1:end-4) '_ssa_Lung' '.wav'],lungSound,fs)


end

%%
save('TempResults\otherMethods.mat')

