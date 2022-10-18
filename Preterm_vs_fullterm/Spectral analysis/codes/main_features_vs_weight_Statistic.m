%%
% This function demonstrate the scatter plot of the extracted featured vs
% the birth weigh, while the calculated Spearman, Pearson and Kendall
% correlation are attached.
% determine the address of feature file at the second line after initialization.
% determine address of the file contatining the weights at the firat line
% after initialization.

% Davood Fattahi, Feb 2020
%% initialization

clear
close all
clc

%% the addresses:
Weight_address='Weights.xlsx';
Feature_addrss = '10secResults2\10secResults-MFCC\MFCC15 features.xls';
% % 
% addrss = '10secResults2\10secResults-VeryHighFreqLin\LogPowerBandFeatures1000-2000Over100-2000.xls';
% addrss =  '10secResults2\10secResults-VeryHighFreqLin\LogPowerBandFeatures1000-2000Over1000-2000.xls';
% addrss = '10secResults2\10secResults-VeryHighFreqLin\PowerBandFeatures1000-2000Over100-2000.xls';
% addrss = '10secResults2\10secResults-VeryHighFreqLin\PowerBandFeatures1000-2000Over1000-2000.xls';
% % 
% addrss = '10secResults2\10secResults-MedFreqLin\LogPowerBandFeatures250-500Over100-1000.xls';
% addrss = '10secResults2\10secResults-MedFreqLin\LogPowerBandFeatures250-500Over250-500.xls';
% addrss = '10secResults2\10secResults-MedFreqLin\PowerBandFeatures250-500Over100-1000.xls';
% addrss = '10secResults2\10secResults-MedFreqLin\PowerBandFeatures250-500Over250-500.xls';
% % 
% addrss = '10secResults2\10secResults-LowFreqLin\LogPowerBandFeatures100-250Over100-250.xls';
% addrss = '10secResults2\10secResults-LowFreqLin\LogPowerBandFeatures100-250Over100-1000.xls';
% addrss = '10secResults2\10secResults-LowFreqLin\PowerBandFeatures100-250Over100-250.xls';
% addrss = '10secResults2\10secResults-LowFreqLin\PowerBandFeatures100-250Over100-1000.xls';
% % 
% addrss = '10secResults2\10secResults-HighFreqLin\LogPowerBandFeatures500-1000Over100-1000.xls';
% addrss = '10secResults2\10secResults-HighFreqLin\LogPowerBandFeatures500-1000Over500-1000.xls';
% addrss = '10secResults2\10secResults-HighFreqLin\PowerBandFeatures500-1000Over100-1000.xls';
% addrss = '10secResults2\10secResults-HighFreqLin\PowerBandFeatures500-1000Over500-1000.xls';
% 

% addrss = '10secResults2\10secResult-MFCCbasedFreqs\LogPowerBandFeatures133-999Over133-1000.xls';
% addrss = '10secResults2\10secResult-MFCCbasedFreqs\PowerBandFeatures133-999Over133-1000.xls';




[Fnum,Ftxt,Fraw]=xlsread( Feature_addrss);
[~,col]=find(Feature_addrss=='\');
FileName=Feature_addrss(col(end)+1:end-4);


[LnumPre,LtxtPre,LrawPre]=xlsread(Weight_address,'Sheet1');
[LnumFul,LtxtFul,LrawFul]=xlsread(Weight_address,'Sheet2');

Fname=Ftxt(2:end,1);
LnamePre=LtxtPre(2:end,1);
LnameFul=LtxtFul(2:end,1);

%% removing extensions form file names
for i=1:size(Fname,1)
    if Fname{i}(end-3:end-2)=='.m'
        Fname{i}(end-3:end)=[];
    end
end
%% Finding the common cases in both labels and features datasheets
[ind,loc]=ismember(LnamePre,Fname);
FnamePre=Fname(loc(loc~=0),:);
FnumPre=Fnum(loc(loc~=0),:);
LnamePre=LnamePre(ind,:);
LnumPre=LnumPre(ind,:);


[ind,loc]=ismember(LnameFul,Fname);
FnameFul=Fname(loc(loc~=0),:);
FnumFul=Fnum(loc(loc~=0),:);
LnameFul=LnameFul(ind,:);
LnumFul=LnumFul(ind,:);

L=[LnumPre;LnumFul];
F=[FnumPre;FnumFul];

CrSp=corr(L,F,'type','Spearman');
CrPr=corr(L,F,'type','Pearson');
CrKn=corr(L,F,'type','Kendall');
PrTrMean=mean(FnumPre);
FlTrMean=mean(FnumFul);
PrTrVar=var(FnumPre);
FlTrVar=var(FnumFul);


for j=1:size(L,2)
    k=0;
    kk=1;
    figure('units','normalized','outerposition',[0 0 1 1])
    for i=1:size(F,2)
        k=k+1;
        if k == 7
            saveas(gcf,[FileName int2str(kk) '.jpeg'])
%             savefig([FileName int2str(kk) '.fig']);
            kk=kk+1;
            figure('units','normalized','outerposition',[0 0 1 1])
            k=k-6;
        end
        subplot(2,3,k)
        scatter(LnumPre(:,j),FnumPre(:,i),'filled','r');
        hold on
        scatter(LnumFul(:,j),FnumFul(:,i),'filled','blue');
        set(gca,'FontSize',6);
        title({['\fontsize{8}Scatter of ' Ftxt{1,1+i} ' and ' LtxtPre{1,j+1}], ...
            ['Spearman=' num2str(CrSp(j,i))], ['Pearson=' num2str(CrPr(j,i))], ...
            ['Kendall=' num2str(CrKn(j,i))], },'FontWeight','normal');
    end
    saveas(gcf,[FileName int2str(kk) '.jpeg'])
%     savefig([FileName int2str(kk) '.fig']);
end
    
