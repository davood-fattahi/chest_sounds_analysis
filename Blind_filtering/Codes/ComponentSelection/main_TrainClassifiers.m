clear
% close all
clc

fs=4000;
load componentSelectionData.mat

%%%% parameter setting for Welch
wind = ones(1,floor(0.06*fs)); % 60 ms
nover = floor(length(wind)*.25); % 0.25 overlap
nfft = 2^(nextpow2(length(wind))-1); % nfft
fcl=50;  %% Low cutoff frequency
fch=2000;    %% High cutoff frequency
fbands=[50 200; 200 400; 400 1000; 1000 1200; 1200 2000];
AHR=2.2;


Features=zeros(size(Cmpnnts,2), size(fbands,1)+3);
[pxx,f] = pwelch(real(Cmpnnts),wind,nover,nfft,fs); % Welch
for i=1:size(Cmpnnts,2)
    bp=bandpower(pxx(:,i),f,[fcl fch],'psd');
    for j=1:size(fbands,1)
        Features(i,j)=bandpower(pxx(:,i),f,[fbands(j,1) fbands(j,2)],'psd')/bp; % Normalized Features
    end
    Peaks=find(PeakDetect(abs(Cmpnnts(:,i)),AHR/fs)); 
    Features(i,j+1)=std(Peaks(2:end)-Peaks(1:end-1));
    Features(i,j+2)=autoCorrPeaksRatio(Cmpnnts(:,i),fs,AHR);
    [~, I]=max(pxx(:,i)); Features(i,j+3)=f(I);
%     Features(i,j+4) = powerbw(pxx,f);
%     Features(i,j+5) = median(pxx); 
end

%%
%%% in the case of feature selection:
% idx = fscmrmr(normalize(Features),Lbls);
% idx=idx(1:6);
%%% else:
idx=1:size(Features,2);

%% model training
rng(1); % fix the random generator seed
Mdl = fitcknn(normalize(Features(:,idx)),Lbls, ...
    'OptimizeHyperparameters','all',...
    'HyperparameterOptimizationOptions',struct('Kfold',10, ...
    'MaxObjectiveEvaluations',100));
save('MdlKnn.mat','Mdl');

rng(1); % fix the random generator seed
Mdl = fitcdiscr(normalize(Features(:,idx)),Lbls, ...
    'OptimizeHyperparameters','auto',...
    'DiscrimType', 'linear', ...
    'HyperparameterOptimizationOptions',struct('Kfold',10, ...
    'MaxObjectiveEvaluations',100));
save('MdlLda.mat','Mdl');

rng(1); % fix the random generator seed
Mdl = fitcecoc(normalize(Features(:,idx)),Lbls, ...
    'OptimizeHyperparameters','all',...
    'HyperparameterOptimizationOptions',struct('Kfold',10, ...
    'MaxObjectiveEvaluations',100));
save('MdlSvm.mat','Mdl');

rng(1); % fix the random generator seed
Mdl = fitcnb(normalize(Features(:,idx)),Lbls, ...
    'OptimizeHyperparameters','all',...
    'HyperparameterOptimizationOptions',struct('Kfold',10, ...
    'MaxObjectiveEvaluations',100));
save('MdlBys.mat','Mdl');

rng(1); % fix the random generator seed
Mdl = fitctree(normalize(Features(:,idx)),Lbls, ...
    'OptimizeHyperparameters','all',...
    'HyperparameterOptimizationOptions',struct('Kfold',10, ...
    'MaxObjectiveEvaluations',100));
save('MdlDecTree.mat','Mdl');


config=struct('wind',wind, 'nover', nover, 'nfft', nfft, 'fcl', fcl, 'fch' ...
    ,fch, 'fbands', fbands, 'AHR', AHR, 'idx', idx );
save('config.mat','config');


