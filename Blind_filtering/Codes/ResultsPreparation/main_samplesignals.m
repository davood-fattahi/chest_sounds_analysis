clear
% close all
clc
% 
% 
% LungInhaleAnn=[0.635 1.000; 1.902 2.294; 2.825 3.155; 4.622 5.081; 5.612 5.970; 6.582 7.208; 7.686 8.317; 9.258 9.424];
% LungExhaleAnn=[0.215 0.535; 1.330 1.683; 2.352 2.734; 3.465 3.795; 5.090 5.339; 5.984 6.424; 7.194 7.614; 8.331 8.722; 9.512 9.923];
% HeartS1S2Ann=[0.215 0.308; 0.421 0.494; 0.625 0.740; 0.913 0.993; 1.109 1.169; ...
%     1.355 1.426; 1.538 1.586; 1.778 1.885; 1.971 2.037; 2.210 2.291; 2.409 2.456; ...
%     2.642 2.714; 2.839 2.895; 3.062 3.144; 3.248 3.326; 3.487 3.551; 3.700 3.749; ...
%     3.918 4.004; 4.111 4.179; 4.339 4.435; 4.555 4.609; 4.774 4.850; 4.981 5.025; ...
%     5.205 5.285; 5.406 5.457; 5.627 5.715; 5.840 5.889; 6.060 6.136; 6.262 6.312; ...
%     6.498 6.571; 6.693 6.753; 6.925 7.002; 7.129 7.185; 7.361 7.433; 7.560 7.611; ...
%     7.786 7.865; 7.991 8.041; 8.220 8.289; 8.426 8.467; 8.660 8.729; 8.851 8.904; ...
%     9.080 9.163; 9.294 9.337; 9.526 9.594; 9.720 9.773; 9.963 10.018];

load samplesignal.mat

LungInhaleInd=false(size(Raw));
for i=1:size(LungInhaleAnn,1)
    LungInhaleInd(LungInhaleAnn(i,1)*fs:LungInhaleAnn(i,2)*fs)=true;
end

LungExhaleInd=false(size(Raw));
for i=1:size(LungExhaleAnn,1)
    LungExhaleInd(LungExhaleAnn(i,1)*fs:LungExhaleAnn(i,2)*fs)=true;
end


HeartS1S2Ind=false(size(Raw));
for i=1:size(HeartS1S2Ann,1)
    HeartS1S2Ind(HeartS1S2Ann(i,1)*fs:(HeartS1S2Ann(i,2)*fs))=true;
end


%%
forgColor=[.25 .35 1 ];
backgColor=[.6 .6 .6];
figure('units','normalized','outerposition',[0 0 .5 1])    
subplot(6,1,1)
tstmp=1/fs:1/fs:size(Raw,2)/fs; tstmp=tstmp-5;
plot(tstmp,Raw,'-','Color', backgColor); xlim([0 5]); ylim([-1 1]);

subplot(6,1,2)
tstmp=1/fs:1/fs:size(FreqFilt_Lung,2)/fs; tstmp=tstmp-5;
plot(tstmp,FreqFilt_Lung,'-','Color', backgColor); xlim([0 5]); ylim([-1 1]);
hold on;
% for i=1:size(LungInhaleAnn,1)
%     indx=LungInhaleAnn(i,1)*fs:LungInhaleAnn(i,2)*fs;
%     plot(tstmp(indx),FreqFilt_Lung(indx),'Color',[.7 .7 .7]);
% end
for i=1:size(LungExhaleAnn,1)
    indx=LungExhaleAnn(i,1)*fs:LungExhaleAnn(i,2)*fs;
    plot(tstmp(indx),FreqFilt_Lung(indx),'Color',forgColor);
end
% for i=1:size(HeartS1S2Ann,1)
%     indx=HeartS1S2Ann(i,1)*fs:HeartS1S2Ann(i,2)*fs;
%     plot(tstmp(indx),FreqFilt_Lung(indx),'Color',[0 0 0]);
% end

subplot(6,1,3)
tstmp=1/fs:1/fs:size(SwtPca_Lung,2)/fs; tstmp=tstmp-5;
plot(tstmp,SwtPca_Lung,'-','Color', backgColor); xlim([0 5]); ylim([-1 1]);
hold on;
% for i=1:size(LungInhaleAnn,1)
%     indx=LungInhaleAnn(i,1)*fs:LungInhaleAnn(i,2)*fs;
%     plot(tstmp(indx),SwtPca_Lung(indx),'Color',forgColor);
% end
for i=1:size(LungExhaleAnn,1)
    indx=LungExhaleAnn(i,1)*fs:LungExhaleAnn(i,2)*fs;
    plot(tstmp(indx),SwtPca_Lung(indx),'Color',forgColor);
end
% for i=1:size(HeartS1S2Ann,1)
%     indx=HeartS1S2Ann(i,1)*fs:HeartS1S2Ann(i,2)*fs;
%     plot(tstmp(indx),SwtPca_Lung(indx),'Color',[0 0 0]);
% end

subplot(6,1,4)
tstmp=1/fs:1/fs:size(CwtPca_Lung,2)/fs; tstmp=tstmp-5;
plot(tstmp,CwtPca_Lung,'-','Color', backgColor); xlim([0 5]); ylim([-1 1]);
hold on;
% for i=1:size(LungInhaleAnn,1)
%     indx=LungInhaleAnn(i,1)*fs:LungInhaleAnn(i,2)*fs;
%     plot(tstmp(indx),CwtPca_Lung(indx),'Color',forgColor);
% end
for i=1:size(LungExhaleAnn,1)
    indx=LungExhaleAnn(i,1)*fs:LungExhaleAnn(i,2)*fs;
    plot(tstmp(indx),CwtPca_Lung(indx),'Color',forgColor);
end
% for i=1:size(HeartS1S2Ann,1)
%     indx=HeartS1S2Ann(i,1)*fs:HeartS1S2Ann(i,2)*fs;
%     plot(tstmp(indx),CwtPca_Lung(indx),'Color',[0 0 0]);
% end

subplot(6,1,5)
tstmp=1/fs:1/fs:size(CwtPica_Lung,2)/fs; tstmp=tstmp-5;
plot(tstmp,CwtPica_Lung,'-','Color', backgColor); xlim([0 5]); ylim([-1 1]);
hold on;
% for i=1:size(LungInhaleAnn,1)
%     indx=LungInhaleAnn(i,1)*fs:LungInhaleAnn(i,2)*fs;
%     plot(tstmp(indx),CwtPica_Lung(indx),'Color',forgColor);
% end
for i=1:size(LungExhaleAnn,1)
    indx=LungExhaleAnn(i,1)*fs:LungExhaleAnn(i,2)*fs;
    plot(tstmp(indx),CwtPica_Lung(indx),'Color',forgColor);
end
% for i=1:size(HeartS1S2Ann,1)
%     indx=HeartS1S2Ann(i,1)*fs:HeartS1S2Ann(i,2)*fs;
%     plot(tstmp(indx),CwtPica_Lung(indx),'Color',[0 0 0]);
% end

subplot(6,1,6)
tstmp=1/fs:1/fs:size(CwtSobi_Lung,2)/fs; tstmp=tstmp-5;
plot(tstmp,CwtSobi_Lung,'-','Color', backgColor); xlim([0 5]); ylim([-1 1]);
hold on;
% for i=1:size(LungInhaleAnn,1)
%     indx=LungInhaleAnn(i,1)*fs:LungInhaleAnn(i,2)*fs;
%     plot(tstmp(indx),CwtSobi_Lung(indx),'Color',forgColor);
% end
for i=1:size(LungExhaleAnn,1)
    indx=LungExhaleAnn(i,1)*fs:LungExhaleAnn(i,2)*fs;
    plot(tstmp(indx),CwtSobi_Lung(indx),'Color',forgColor);
end
% for i=1:size(HeartS1S2Ann,1)
%     indx=HeartS1S2Ann(i,1)*fs:HeartS1S2Ann(i,2)*fs;
%     plot(tstmp(indx),CwtSobi_Lung(indx),'Color',[0 0 0]);
% end
savefig('LungSignals.fig');
saveas(gcf,'LungSignals.eps','epsc');

%%
forgColor=[1 .3 .3 ];
backgColor=[.6 .6 .6];
figure('units','normalized','outerposition',[0 0 .5 1])    
subplot(6,1,1)
tstmp=1/fs:1/fs:size(Raw,2)/fs; tstmp=tstmp-5;
plot(tstmp,Raw,'-','Color', backgColor ); xlim([0 5]); ylim([-1 1]);

subplot(6,1,2)
tstmp=1/fs:1/fs:size(FreqFilt_Heart,2)/fs; tstmp=tstmp-5;
plot(tstmp,FreqFilt_Heart,'-','Color', backgColor); xlim([0 5]); ylim([-1 1]);
hold on;
for i=1:size(HeartS1S2Ann,1)
    indx=HeartS1S2Ann(i,1)*fs:HeartS1S2Ann(i,2)*fs;
    plot(tstmp(indx),FreqFilt_Heart(indx),'Color',forgColor);
end

subplot(6,1,3)
tstmp=1/fs:1/fs:size(SwtPca_Heart,2)/fs; tstmp=tstmp-5;
plot(tstmp,SwtPca_Heart,'-','Color', backgColor); xlim([0 5]); ylim([-1 1]);
hold on;
for i=1:size(HeartS1S2Ann,1)
    indx=HeartS1S2Ann(i,1)*fs:HeartS1S2Ann(i,2)*fs;
    plot(tstmp(indx),SwtPca_Heart(indx),'Color',forgColor);
end

subplot(6,1,4)
tstmp=1/fs:1/fs:size(CwtPca_Heart,2)/fs; tstmp=tstmp-5;
plot(tstmp,CwtPca_Heart,'-','Color', backgColor); xlim([0 5]); ylim([-1 1]);
hold on;
for i=1:size(HeartS1S2Ann,1)
    indx=HeartS1S2Ann(i,1)*fs:HeartS1S2Ann(i,2)*fs;
    plot(tstmp(indx),CwtPca_Heart(indx),'Color',forgColor);
end

subplot(6,1,5)
tstmp=1/fs:1/fs:size(CwtPica_Heart,2)/fs; tstmp=tstmp-5;
plot(tstmp,CwtPica_Heart,'-','Color', backgColor); xlim([0 5]); ylim([-1 1]);
hold on;
for i=1:size(HeartS1S2Ann,1)
    indx=HeartS1S2Ann(i,1)*fs:HeartS1S2Ann(i,2)*fs;
    plot(tstmp(indx),CwtPica_Heart(indx),'Color',forgColor);
end

subplot(6,1,6)
tstmp=1/fs:1/fs:size(CwtSobi_Heart,2)/fs; tstmp=tstmp-5;
plot(tstmp,CwtSobi_Heart,'-','Color', backgColor); xlim([0 5]); ylim([-1 1]);
hold on;
for i=1:size(HeartS1S2Ann,1)
    indx=HeartS1S2Ann(i,1)*fs:HeartS1S2Ann(i,2)*fs;
    plot(tstmp(indx),CwtSobi_Heart(indx),'Color',forgColor);
end
savefig('HeartSignals.fig');
saveas(gcf,'HeartSignals.eps','epsc');


