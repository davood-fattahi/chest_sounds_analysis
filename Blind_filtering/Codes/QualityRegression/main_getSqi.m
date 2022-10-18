clear
% close all

all=readcell('rawHeartAllFeatures.xlsx');
ftrs=sortrows((all(:,2:end))')'; 
ftrs=cell2mat(ftrs(2:end,:));
load('heartAllFeatCVMdl.mat');
rawHeartSqi = predict(heartCVMdl,normalize(ftrs(:,heartFeatSelectIdx)')');
names=all(2:end,1); [~,I]=sort(names); 
rawHeartSqi=rawHeartSqi(I);


all=readcell('freqFixHeartAllFeatures.xlsx');
ftrs=sortrows((all(:,2:end))')'; 
ftrs=cell2mat(ftrs(2:end,:));
load('heartAllFeatCVMdl.mat');
freqFixHeartSqi = predict(heartCVMdl,normalize(ftrs(:,heartFeatSelectIdx)')');
names=all(2:end,1);
names = cellfun(@(x) [x(1:regexp(x,'_Freq')-1) '.wav'], names,'UniformOutput',false);
[~,I]=sort(names); 
freqFixHeartSqi=freqFixHeartSqi(I);

all=readcell('swtPcaHeartAllFeatures.xlsx');
ftrs=sortrows((all(:,2:end))')'; 
ftrs=cell2mat(ftrs(2:end,:));
load('heartAllFeatCVMdl.mat');
swtPcaHeartSqi = predict(heartCVMdl,normalize(ftrs(:,heartFeatSelectIdx)')');
names=all(2:end,1);
names = cellfun(@(x) [x(1:regexp(x,'_Sw')-1) '.wav'], names,'UniformOutput',false);
[~,I]=sort(names); 
swtPcaHeartSqi=swtPcaHeartSqi(I);


all=readcell('cwtPcaHeartAllFeatures.xlsx');
ftrs=sortrows((all(:,2:end))')'; 
ftrs=cell2mat(ftrs(2:end,:));
load('heartAllFeatCVMdl.mat');
cwtPcaHeartSqi = predict(heartCVMdl,normalize(ftrs(:,heartFeatSelectIdx)')');
names=all(2:end,1);
names = cellfun(@(x) [x(1:regexp(x,'_Cw')-1) '.wav'], names,'UniformOutput',false);
[~,I]=sort(names); 
cwtPcaHeartSqi=cwtPcaHeartSqi(I);


all=readcell('cwtPicaHeartAllFeatures.xlsx');
ftrs=sortrows((all(:,2:end))')'; 
ftrs=cell2mat(ftrs(2:end,:));
load('heartAllFeatCVMdl.mat');
cwtPicaHeartSqi = predict(heartCVMdl,normalize(ftrs(:,heartFeatSelectIdx)')');
names=all(2:end,1);
names = cellfun(@(x) [x(1:regexp(x,'_Cw')-1) '.wav'], names,'UniformOutput',false);
[~,I]=sort(names); 
cwtPicaHeartSqi=cwtPicaHeartSqi(I);


all=readcell('cwtSobiHeartAllFeatures.xlsx');
ftrs=sortrows((all(:,2:end))')'; 
ftrs=cell2mat(ftrs(2:end,:));
load('heartAllFeatCVMdl.mat');
cwtSobiHeartSqi = predict(heartCVMdl,normalize(ftrs(:,heartFeatSelectIdx)')');
names=all(2:end,1);
names = cellfun(@(x) [x(1:regexp(x,'_Cw')-1) '.wav'], names,'UniformOutput',false);
[~,I]=sort(names); 
cwtSobiHeartSqi=cwtSobiHeartSqi(I);


HeartSqiAllF=table(rawHeartSqi, freqFixHeartSqi, swtPcaHeartSqi, cwtPcaHeartSqi, ...
    cwtPicaHeartSqi, cwtSobiHeartSqi);
writetable(HeartSqiAllF, 'HeartSqiAllF.xlsx');


%%

all=readcell('rawLungAllFeatures.xlsx');
ftrs=sortrows((all(:,2:end))')'; 
ftrs=cell2mat(ftrs(2:end,:));
load('lungAllFeatCVMdl.mat');
rawLungSqi = predict(lungCVMdl,normalize(ftrs(:,lungFeatSelectIdx)')');
names=all(2:end,1); [~,I]=sort(names); 
rawLungSqi=rawLungSqi(I);


all=readcell('freqFixLungAllFeatures.xlsx');
ftrs=sortrows((all(:,2:end))')'; 
ftrs=cell2mat(ftrs(2:end,:));
load('lungAllFeatCVMdl.mat');
freqFixLungSqi = predict(lungCVMdl,normalize(ftrs(:,lungFeatSelectIdx)')');
names=all(2:end,1);
names = cellfun(@(x) [x(1:regexp(x,'_Freq')-1) '.wav'], names,'UniformOutput',false);
[~,I]=sort(names); 
freqFixLungSqi=freqFixLungSqi(I);

all=readcell('swtPcaLungAllFeatures.xlsx');
ftrs=sortrows((all(:,2:end))')'; 
ftrs=cell2mat(ftrs(2:end,:));
load('lungAllFeatCVMdl.mat');
swtPcaLungSqi = predict(lungCVMdl,normalize(ftrs(:,lungFeatSelectIdx)')');
names=all(2:end,1);
names = cellfun(@(x) [x(1:regexp(x,'_Sw')-1) '.wav'], names,'UniformOutput',false);
[~,I]=sort(names); 
swtPcaLungSqi=swtPcaLungSqi(I);


all=readcell('cwtPcaLungAllFeatures.xlsx');
ftrs=sortrows((all(:,2:end))')'; 
ftrs=cell2mat(ftrs(2:end,:));
load('lungAllFeatCVMdl.mat');
cwtPcaLungSqi = predict(lungCVMdl,normalize(ftrs(:,lungFeatSelectIdx)')');
names=all(2:end,1);
names = cellfun(@(x) [x(1:regexp(x,'_Cw')-1) '.wav'], names,'UniformOutput',false);
[~,I]=sort(names); 
cwtPcaLungSqi=cwtPcaLungSqi(I);


all=readcell('cwtPicaLungAllFeatures.xlsx');
ftrs=sortrows((all(:,2:end))')'; 
ftrs=cell2mat(ftrs(2:end,:));
load('lungAllFeatCVMdl.mat');
cwtPicaLungSqi = predict(lungCVMdl,normalize(ftrs(:,lungFeatSelectIdx)')');
names=all(2:end,1);
names = cellfun(@(x) [x(1:regexp(x,'_Cw')-1) '.wav'], names,'UniformOutput',false);
[~,I]=sort(names); 
cwtPicaLungSqi=cwtPicaLungSqi(I);


all=readcell('cwtSobiLungAllFeatures.xlsx');
ftrs=sortrows((all(:,2:end))')'; 
ftrs=cell2mat(ftrs(2:end,:));
load('lungAllFeatCVMdl.mat');
cwtSobiLungSqi = predict(lungCVMdl,normalize(ftrs(:,lungFeatSelectIdx)')');
names=all(2:end,1);
names = cellfun(@(x) [x(1:regexp(x,'_Cw')-1) '.wav'], names,'UniformOutput',false);
[~,I]=sort(names); 
cwtSobiLungSqi=cwtSobiLungSqi(I);


LungSqiAllF=table(rawLungSqi, freqFixLungSqi, swtPcaLungSqi, cwtPcaLungSqi, ...
    cwtPicaLungSqi, cwtSobiLungSqi);
writetable(LungSqiAllF, 'LungSqiAllF.xlsx');


%% Other methods ...

all=readcell('emdHeartAllFeatures.xlsx');
ftrs=sortrows((all(:,2:end))')'; 
ftrs=cell2mat(ftrs(2:end,:));
load('heartAllFeatCVMdl.mat');
emdHeartSqi = predict(heartCVMdl,normalize(ftrs(:,heartFeatSelectIdx)')');
names=all(2:end,1);
names = cellfun(@(x) [x(1:regexp(x,'_emd')-1) '.wav'], names,'UniformOutput',false);
[~,I]=sort(names); 
emdHeartSqi=emdHeartSqi(I);

all=readcell('emdLungAllFeatures.xlsx');
ftrs=sortrows((all(:,2:end))')'; 
ftrs=cell2mat(ftrs(2:end,:));
load('lungAllFeatCVMdl.mat');
emdLungSqi = predict(lungCVMdl,normalize(ftrs(:,lungFeatSelectIdx)')');
names=all(2:end,1);
names = cellfun(@(x) [x(1:regexp(x,'_emd')-1) '.wav'], names,'UniformOutput',false);
[~,I]=sort(names); 
emdLungSqi=emdLungSqi(I);


%%%%%%%%%%%%%%%%%%%%%

all=readcell('eemdHeartAllFeatures.xlsx');
ftrs=sortrows((all(:,2:end))')'; 
ftrs=cell2mat(ftrs(2:end,:));
load('heartAllFeatCVMdl.mat');
eemdHeartSqi = predict(heartCVMdl,normalize(ftrs(:,heartFeatSelectIdx)')');
names=all(2:end,1);
names = cellfun(@(x) [x(1:regexp(x,'_eemd')-1) '.wav'], names,'UniformOutput',false);
[~,I]=sort(names); 
eemdHeartSqi=eemdHeartSqi(I);

all=readcell('eemdLungAllFeatures.xlsx');
ftrs=sortrows((all(:,2:end))')'; 
ftrs=cell2mat(ftrs(2:end,:));
load('lungAllFeatCVMdl.mat');
eemdLungSqi = predict(lungCVMdl,normalize(ftrs(:,lungFeatSelectIdx)')');
names=all(2:end,1);
names = cellfun(@(x) [x(1:regexp(x,'_eemd')-1) '.wav'], names,'UniformOutput',false);
[~,I]=sort(names); 
eemdLungSqi=eemdLungSqi(I);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

all=readcell('ceemdHeartAllFeatures.xlsx');
ftrs=sortrows((all(:,2:end))')'; 
ftrs=cell2mat(ftrs(2:end,:));
load('heartAllFeatCVMdl.mat');
ceemdHeartSqi = predict(heartCVMdl,normalize(ftrs(:,heartFeatSelectIdx)')');
names=all(2:end,1);
names = cellfun(@(x) [x(1:regexp(x,'_ceemd')-1) '.wav'], names,'UniformOutput',false);
[~,I]=sort(names); 
ceemdHeartSqi=ceemdHeartSqi(I);

all=readcell('ceemdLungAllFeatures.xlsx');
ftrs=sortrows((all(:,2:end))')'; 
ftrs=cell2mat(ftrs(2:end,:));
load('lungAllFeatCVMdl.mat');
ceemdLungSqi = predict(lungCVMdl,normalize(ftrs(:,lungFeatSelectIdx)')');
names=all(2:end,1);
names = cellfun(@(x) [x(1:regexp(x,'_ceemd')-1) '.wav'], names,'UniformOutput',false);
[~,I]=sort(names); 
ceemdLungSqi=emdLungSqi(I);



%%%%%%%%%%%%%%%%%%%%%%%%%%

all=readcell('ssaHeartAllFeatures.xlsx');
ftrs=sortrows((all(:,2:end))')'; 
ftrs=cell2mat(ftrs(2:end,:));
load('heartAllFeatCVMdl.mat');
ssaHeartSqi = predict(heartCVMdl,normalize(ftrs(:,heartFeatSelectIdx)')');
names=all(2:end,1);
names = cellfun(@(x) [x(1:regexp(x,'_ssa')-1) '.wav'], names,'UniformOutput',false);
[~,I]=sort(names); 
ssaHeartSqi=ssaHeartSqi(I);
ssaHeartSqi(ssaHeartSqi>5)=5; 
ssaHeartSqi(isnan(ssaHeartSqi))=1;

all=readcell('ssaLungAllFeatures.xlsx');
ftrs=sortrows((all(:,2:end))')'; 
ftrs=cell2mat(ftrs(2:end,:));
load('lungAllFeatCVMdl.mat');
ssaLungSqi = predict(lungCVMdl,normalize(ftrs(:,lungFeatSelectIdx)')');
names=all(2:end,1);
names = cellfun(@(x) [x(1:regexp(x,'_ssa')-1) '.wav'], names,'UniformOutput',false);
[~,I]=sort(names); 
ssaLungSqi=ssaLungSqi(I);
ssaLungSqi(ssaLungSqi>5)=5; 
ssaLungSqi(isnan(ssaLungSqi))=1;
%%%%%%%%%%%%%%%%%%%%%%%%%

HeartSqiAllFeatOtherMethods=table(emdHeartSqi, eemdHeartSqi, ceemdHeartSqi, ssaHeartSqi);
writetable(HeartSqiAllFeatOtherMethods, 'HeartSqiAllFeatOtherMethods.xlsx');

LungSqiAllFeatOtherMethods = table(emdLungSqi, eemdLungSqi, ceemdLungSqi, ssaLungSqi);
writetable(LungSqiAllFeatOtherMethods, 'LungSqiAllFeatOtherMethods.xlsx');





