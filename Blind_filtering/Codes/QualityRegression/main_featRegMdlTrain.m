
clear
% close all


%% get features
heartAllFtrs=readcell('trainHeartAllFeatures.xlsx');
lungAllFtrs=readcell('trainLungAllFeatures.xlsx');

heartFtrsIdx=[];
lungFtrsIdx=[];


heartAllFtrs=sortrows((heartAllFtrs(:,2:end))')';
lungAllFtrs=sortrows((lungAllFtrs(:,2:end))')';
heartAllFtrs=cell2mat(heartAllFtrs(2:end,:));
lungAllFtrs=cell2mat(lungAllFtrs(2:end,:));

%% get annotations
heartQul5An=readcell('heartQul5Annots.xlsx');
heartQul5An=cell2mat(heartQul5An(:,2:end));
lungQul5An=readcell('lungQul5Annots.xlsx');
lungQul5An=cell2mat(lungQul5An(:,2:end));

%% regression training
close all

rng(1)
heartFeatSelectIdx =  fscmrmr(normalize(heartAllFtrs')',mean(heartQul5An,2));
heartFeatSelectIdx=heartFeatSelectIdx(1:40);
heartCVMdl = fitrlinear(normalize(heartAllFtrs(:,heartFeatSelectIdx)')',mean(heartQul5An,2), ...
    'OptimizeHyperparameters', 'all', 'HyperparameterOptimizationOptions', ...
    struct('AcquisitionFunctionName','expected-improvement-plus', ...
    'MaxObjectiveEvaluations',100));
% save('heartAllFeatCVMdl.mat','heartCVMdl', 'heartFeatSelectIdx')

rng(1)
lungFeatSelectIdx = fscmrmr(normalize(lungAllFtrs')',mean(lungQul5An,2));
lungFeatSelectIdx=lungFeatSelectIdx(1:40);
lungCVMdl = fitrlinear(normalize(lungAllFtrs(:,lungFeatSelectIdx)')',mean(lungQul5An,2), ...
    'OptimizeHyperparameters', 'all', 'HyperparameterOptimizationOptions', ...
    struct('AcquisitionFunctionName','expected-improvement-plus', ...
    'MaxObjectiveEvaluations',100));

% save('lungAllFeatCVMdl.mat','lungCVMdl', 'lungFeatSelectIdx')

%% other checked models
% heartFeatSelectIdx =  fscchi2(normalize(heartAllFtrs')',mean(heartQul5An,2));
% heartFeatSelectIdx = fsrftest(normalize(heartAllFtrs')',mean(heartQul5An,2));
% [heartFeatSelectIdx, weights] = relieff(normalize(heartAllFtrs')',mean(heartQul5An,2), 3);
% heartFeatSelectIdx =  fscmrmr(normalize(heartAllFtrs')',mean(heartQul5An,2));


% heartFeatSelectIdx=heartFeatSelectIdx(1:40);
% heartCVMdl = fitrgp(normalize(heartAllFtrs(:,heartFeatSelectIdx)')',mean(heartQul5An,2), ...
%     'OptimizeHyperparameters','all' , 'HyperparameterOptimizationOptions', ...
%     struct('AcquisitionFunctionName','expected-improvement-plus', ...
%     'MaxObjectiveEvaluations',20));






