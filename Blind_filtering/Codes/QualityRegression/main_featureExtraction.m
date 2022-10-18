clc 
clear 
% close all

%% train set features
% load annotations
annotations.Heart=readcell("Summary of All Annotations (BFF paper).xlsx", 'Sheet',"Heart");
% load directory information of the files
files.Heart=dir("Annotation\*\Raw Files\Heart");
% for each annoatted heart record:
filesName=annotations.Heart(2:end,strcmp(annotations.Heart(1,:),'File'));
DA={files.Heart.folder}'; DA=DA(matches({files.Heart.name}',filesName));
trainHeartAllFeatures=collectAllFeatures(DA,filesName);
writetable(trainHeartAllFeatures, 'trainHeartAllFeatures.xlsx','WriteRowNames',true,'WriteVariableNames',true);

% get annotations
heartQul5Annots=sort(string((annotations.Heart(2:end,7:13))),2);
heartQul5Annots=[filesName num2cell(str2double(heartQul5Annots(:,1:5)))];
writecell(heartQul5Annots, 'heartQul5Annots.xlsx');

% load annotations
annotations.Lung=readcell("Summary of All Annotations (BFF paper).xlsx", 'Sheet',"Lung");
% load directory information of the files
files.Lung=dir("Annotation\*\Raw Files\Lung");
% for each annoatted lung record:
filesName=annotations.Lung(2:end,strcmp(annotations.Lung(1,:),'File'));
DA={files.Lung.folder}'; DA=DA(matches({files.Lung.name}',filesName));
trainLungAllFeatures=collectAllFeatures(DA,filesName);
writetable(trainLungAllFeatures, 'trainLungAllFeatures.xlsx','WriteRowNames',true,'WriteVariableNames',true);

% get annotations
lungQul5Annots=sort(string((annotations.Lung(2:end,7:13))),2);
lungQul5Annots=[filesName num2cell(str2double(lungQul5Annots(:,1:5)))];
writecell(lungQul5Annots, 'lungQul5Annots.xlsx');

% load tempsave.mat
% save('tempsave.mat');


%% raw records features --------------------------------------------------
DA='..\..\ChestRecords\CurrentDataset\';
dataSheet=readcell([DA '00dataSheet.xlsx']);
filesName=dataSheet([false; cell2mat(dataSheet(2:end,5))==1],1);
rawHeartAllFeatures=collectAllFeatures(DA,filesName);
rawLungAllFeatures=rawHeartAllFeatures;
writetable(rawHeartAllFeatures, 'rawHeartAllFeatures.xlsx','WriteRowNames',true,'WriteVariableNames',true);
writetable(rawLungAllFeatures, 'rawLungAllFeatures.xlsx','WriteRowNames',true,'WriteVariableNames',true);
% save('tempsave.mat');

%% Fix freq features -----------------------------------------------------
DA='..\FilteringMethods\FrequencyFiltering\Results\FreqFixHeart\';
files=dir([DA '\*.wav']);
freqFixHeartAllFeatures=collectAllFeatures(DA,{files.name}');
DA='..\FilteringMethods\FrequencyFiltering\Results\FreqFixLung\';
files=dir([DA '\*.wav']); 
freqFixLungAllFeatures=collectAllFeatures(DA,{files.name}');
writetable(freqFixHeartAllFeatures, 'freqFixHeartAllFeatures.xlsx','WriteRowNames',true,'WriteVariableNames',true);
writetable(freqFixLungAllFeatures, 'freqFixLungAllFeatures.xlsx','WriteRowNames',true,'WriteVariableNames',true);
% save('tempsave.mat');

%% SwtPca features -----------------------------------------------
DA='..\FilteringMethods\swt-pca\TempResults\SwtPcaHeart\';
files=dir([DA '\*.wav']); 
[swtPcaHeartAllFeatures]=collectAllFeatures(DA,{files.name}');
DA='..\FilteringMethods\swt-pca\TempResults\SwtPcaLung\';
files=dir([DA '\*.wav']); 
[swtPcaLungAllFeatures]=collectAllFeatures(DA,{files.name}');
writetable(swtPcaHeartAllFeatures, 'swtPcaHeartAllFeatures.xlsx','WriteRowNames',true,'WriteVariableNames',true);
writetable(swtPcaLungAllFeatures, 'swtPcaLungAllFeatures.xlsx','WriteRowNames',true,'WriteVariableNames',true);
% save('tempsave.mat');

%% CwtPca features -----------------------------------------------
DA='..\FilteringMethods\cwt-pca\TempResults\CwtPcaHeart\';
files=dir([DA '\*.wav']); 
[cwtPcaHeartAllFeatures]=collectAllFeatures(DA,{files.name}');
DA='..\FilteringMethods\cwt-pca\TempResults\CwtPcaLung\';
files=dir([DA '\*.wav']); 
[cwtPcaLungAllFeatures]=collectAllFeatures(DA,{files.name}');
writetable(cwtPcaHeartAllFeatures, 'cwtPcaHeartAllFeatures.xlsx','WriteRowNames',true,'WriteVariableNames',true);
writetable(cwtPcaLungAllFeatures, 'cwtPcaLungAllFeatures.xlsx','WriteRowNames',true,'WriteVariableNames',true);
% save('tempsave.mat');

%% CwtPica features -----------------------------------------------
DA='..\FilteringMethods\cwt-pica\TempResults\CwtPicaHeart\';
files=dir([DA '\*.wav']); 
[cwtPicaHeartAllFeatures]=collectAllFeatures(DA,{files.name}');
DA='..\FilteringMethods\cwt-pica\TempResults\CwtPicaLung\';
files=dir([DA '\*.wav']); 
[cwtPicaLungAllFeatures]=collectAllFeatures(DA,{files.name}');
writetable(cwtPicaHeartAllFeatures, 'cwtPicaHeartAllFeatures.xlsx','WriteRowNames',true,'WriteVariableNames',true);
writetable(cwtPicaLungAllFeatures, 'cwtPicaLungAllFeatures.xlsx','WriteRowNames',true,'WriteVariableNames',true);
% save('tempsave.mat');

%% CwtSobi features -----------------------------------------------
DA='..\FilteringMethods\cwt-sobi\TempResults\CwtSobiHeart\';
files=dir([DA '\*.wav']); 
[cwtSobiHeartAllFeatures]=collectAllFeatures(DA,{files.name}');
DA='..\FilteringMethods\cwt-sobi\TempResults\CwtSobiLung\';
files=dir([DA '\*.wav']); 
[cwtSobiLungAllFeatures]=collectAllFeatures(DA,{files.name}');
writetable(cwtSobiHeartAllFeatures, 'cwtSobiHeartAllFeatures.xlsx','WriteRowNames',true,'WriteVariableNames',true);
writetable(cwtSobiLungAllFeatures, 'cwtSobiLungAllFeatures.xlsx','WriteRowNames',true,'WriteVariableNames',true);
% save('tempsave.mat');

%% EMD features -----------------------------------------------
DA='..\FilteringMethods\otherMethods\TempResults\emdHeart\';
files=dir([DA '\*.wav']); 
[emdHeartAllFeatures]=collectAllFeatures(DA,{files.name}');
DA='..\FilteringMethods\otherMethods\TempResults\emdLung\';
files=dir([DA '\*.wav']); 
[emdLungAllFeatures]=collectAllFeatures(DA,{files.name}');
writetable(emdHeartAllFeatures, 'emdHeartAllFeatures.xlsx','WriteRowNames',true,'WriteVariableNames',true);
writetable(emdLungAllFeatures, 'emdLungAllFeatures.xlsx','WriteRowNames',true,'WriteVariableNames',true);
% save('tempsave.mat');

%% EEMD features -----------------------------------------------
DA='..\FilteringMethods\otherMethods\TempResults\eemdHeart\';
files=dir([DA '\*.wav']); 
[eemdHeartAllFeatures]=collectAllFeatures(DA,{files.name}');
DA='..\FilteringMethods\otherMethods\TempResults\eemdLung\';
files=dir([DA '\*.wav']); 
[eemdLungAllFeatures]=collectAllFeatures(DA,{files.name}');
writetable(eemdHeartAllFeatures, 'eemdHeartAllFeatures.xlsx','WriteRowNames',true,'WriteVariableNames',true);
writetable(eemdLungAllFeatures, 'eemdLungAllFeatures.xlsx','WriteRowNames',true,'WriteVariableNames',true);
% save('tempsave.mat');

%% CEEMD features -----------------------------------------------
DA='..\FilteringMethods\otherMethods\TempResults\ceemdHeart\';
files=dir([DA '\*.wav']); 
[ceemdHeartAllFeatures]=collectAllFeatures(DA,{files.name}');
DA='..\FilteringMethods\otherMethods\TempResults\ceemdLung\';
files=dir([DA '\*.wav']); 
[ceemdLungAllFeatures]=collectAllFeatures(DA,{files.name}');
writetable(ceemdHeartAllFeatures, 'ceemdHeartAllFeatures.xlsx','WriteRowNames',true,'WriteVariableNames',true);
writetable(ceemdLungAllFeatures, 'ceemdLungAllFeatures.xlsx','WriteRowNames',true,'WriteVariableNames',true);
% save('tempsave.mat');


%% ssa features -----------------------------------------------
DA='..\FilteringMethods\otherMethods\TempResults\ssaHeart\';
files=dir([DA '\*.wav']); 
[ssaHeartAllFeatures]=collectAllFeatures(DA,{files.name}');
DA='..\FilteringMethods\otherMethods\TempResults\ssaLung\';
files=dir([DA '\*.wav']); 
[ssaLungAllFeatures]=collectAllFeatures(DA,{files.name}');
writetable(ssaHeartAllFeatures, 'ssaHeartAllFeatures.xlsx','WriteRowNames',true,'WriteVariableNames',true);
writetable(ssaLungAllFeatures, 'ssaLungAllFeatures.xlsx','WriteRowNames',true,'WriteVariableNames',true);
% save('tempsave.mat');



