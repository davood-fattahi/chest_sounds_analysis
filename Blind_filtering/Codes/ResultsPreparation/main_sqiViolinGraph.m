clear
close all
clc

%%% Load Quality Analysis
ourMethodsSqi_Heart=readtable('..\QualityRegression\HeartSqiAllF.xlsx');
ourMethodsSqi_Lung =readtable('..\QualityRegression\LungSqiAllF.xlsx');
otherMethodsSqi_Heart= readtable('..\QualityRegression\HeartSqiAllFeatOtherMethods.xlsx');
otherMethodsSqi_Lung= readtable('..\QualityRegression\LungSqiAllFeatOtherMethods.xlsx');
nmfMethodsSqi_Heart= readtable('..\FilteringMethods\otherMethods\HeartSqiNmfMethods.xlsx');
nmfMethodsSqi_Lung= readtable('..\FilteringMethods\otherMethods\LungSqiNmfMethods.xlsx');
QA_Heart=[ourMethodsSqi_Heart otherMethodsSqi_Heart nmfMethodsSqi_Heart];
QA_Lung=[ourMethodsSqi_Lung otherMethodsSqi_Lung nmfMethodsSqi_Lung];

methodsHeart=["rawHeartSqi"; "freqFixHeartSqi"; "swtPcaHeartSqi"; "cwtPcaHeartSqi"; "cwtPicaHeartSqi"; "cwtSobiHeartSqi"; "ssaHeartSqi"; ...
    "emdHeartSqi"; "eemdHeartSqi"; "ceemdHeartSqi"; "nmfc_STFT_nmf_kl_Filtering"; ...
    "nmfc_STFT_nmf_general_Synthesis"; "nmfc_STFT_nmf_sparse_Filtering" ];

methodsLung=["rawLungSqi"; "freqFixLungSqi"; "swtPcaLungSqi"; "cwtPcaLungSqi"; "cwtPicaLungSqi"; "cwtSobiLungSqi"; "ssaLungSqi"; ...
    "emdLungSqi"; "eemdLungSqi"; "ceemdLungSqi"; "nmfc_STFT_nmf_kl_Filtering"; ...
    "nmfc_STFT_nmf_general_Synthesis"; "nmfc_STFT_nmf_sparse_Filtering" ];

QA_Heart = QA_Heart(:,methodsHeart);
QA_Lung = QA_Lung(:,methodsLung);

QA_Heart=table2array(QA_Heart);
deltaQ_Heart=(QA_Heart(:,2:end))-(QA_Heart(:,1));

QA_Lung=table2array(QA_Lung);
deltaQ_Lung=(QA_Lung(:,2:end))-(QA_Lung(:,1));


MH=mean(deltaQ_Heart);
ML=mean(deltaQ_Lung);

SH=std(deltaQ_Heart);
SL= std(deltaQ_Lung);


%%
figure('Units','normalized','Position',[0.2 0.2 .5 .5]);
boxplot(deltaQ_Heart,'OutlierSize',0.01); hold on;
[h,L,~,~]=violin(deltaQ_Heart,'x',[1 2 3 4 5 6 7 8 9 10 11 12],'facecolor',[.5 .5 .5], ...
    'facealpha', .4, 'edgecolor','none',...
    'bw',0.2,'mc','','medc','','plotlegend',[]);
errorbar((1:12),MH,SH,'.','MarkerSize',10,...
    'Color','k','LineWidth',1); xticks(1:12)

set(gca,'xticklabel',{'Freq. Filt.','SWT-PCA', 'CWT-PCA', 'CWT-\piCA', 'CWT-SOBI', ...
    'SSA', 'EMD', 'EEMD', 'CEEMD', 'NMFC-KL', 'NMFC-L2', 'NMFC-Sparse' }, ...
    'TickLabelInterpreter', 'tex', 'FontSize',10)

xtickangle(60)
ylabel('\Delta quality'); 
set(L,'FontSize',8)
set(h,'LineWidth',.3, 'LineStyle', '--')
xlim([0 13]); ylim([-1 3.5]);
grid on
grid minor

saveas(gcf,'HeartSoundQualityImprovementViolin.fig')
saveas(gcf,'HeartSoundQualityImprovementViolin.eps','epsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure('Units','normalized','Position',[0.2 0.2 .5 .5]);
boxplot(deltaQ_Lung,'OutlierSize',0.01); hold on;
[h,L,~,~]=violin(deltaQ_Lung,'x',[1 2 3 4 5 6 7 8 9 10 11 12],'facecolor',[.5 .5 .5], ...
    'facealpha', .4, 'edgecolor','none',...
    'bw',0.2,'mc','','medc','','plotlegend',[]);
errorbar((1:12),ML,SL,'.','MarkerSize',10,...
    'Color','k','LineWidth',1); xticks(1:12)

set(gca,'xticklabel',{'Freq. Filt.','SWT-PCA', 'CWT-PCA', 'CWT-\piCA', 'CWT-SOBI', ...
    'SSA', 'EMD', 'EEMD', 'CEEMD', 'NMFC-KL', 'NMFC-L2', 'NMFC-Sparse' }, ...
    'TickLabelInterpreter', 'tex', 'FontSize',10)

xtickangle(60)
ylabel('\Delta quality'); 
set(L,'FontSize',8)
set(h,'LineWidth',.3, 'LineStyle', '--')
xlim([0 13]); ylim([-2 3.5]);
grid on
grid minor


saveas(gcf,'LungSoundQualityImprovementViolin.fig')
saveas(gcf,'LungSoundQualityImprovementViolin.eps','epsc')





