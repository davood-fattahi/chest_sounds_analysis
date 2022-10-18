%%%%%%%%%%%%%%%%%%%
%% Initializing
clear
% close all
clc
% delete 'TempResults\CwtPicaHeart\*.wav'
% delete 'TempResults\CwtPicaLung\*.wav'

% Defining address of the folder containing the audio records
DA='..\..\..\chestRecords\currentDataset\';
dataSheet=readcell([DA '00dataSheet.xlsx']);
filesName=dataSheet([false; cell2mat(dataSheet(2:end,5))==1],1);

%% parameter setting (default values)
fl=50; % lower frequency limit
fh=900; % higher frequency limit
vpo=4; % voice per octave
wname='amor'; % wavelet name
AHR=2.2; % average heart rate (Hz);
fpd=75; % approximate frequency for s1 peak detection

%% load the componenet selection model
load config.mat
load 'MdlKnn.mat'

%% do for each of the records
for i= 1 :size(filesName,1)
    %%% Display the progress percentage 
    clc;disp([num2str(i) '/' num2str(size(filesName,1)) ', ' num2str(floor(100*i/size(filesName,1))) ' %']) 

    %%% Loading the signal
    [data,fs]=audioread([DA  '\' filesName{i}]);    %% read the record
    
    %%% Preprocessing - Downsampling to 4000 Hz
    data=downsample(data,fs./4000); fs=4000;
    ff=0:fs/2/length(data):fs/2;
    
    %%% setting length to 10 sec
    data=data(1:10*fs);
    
    %%% Preprocessing - normalizing
    data=normalize(data);  
    
    
    %% cwtpica
    %%% s1 peak detection
    %%% filterbank generation
    fb = cwtfilterbank('SignalLength',length(data),'Wavelet',wname, ...
    'VoicesPerOctave',vpo, 'SamplingFrequency',fs, ...
    'FrequencyLimits',[fl fh]); 

    %%% cwt coefficients (complex value) 
    [Cx, f]=cwt(data,'FilterBank',fb);   
 
    %%% s1 peak detection
    k=find((f<fpd),1);
    Peaks1=PeakDetect(abs(Cx(k,:)),AHR/fs);
    
    
%     figure; plot(real(Cx(end,:))); hold on; plot(Peaks1);   
    [S, Fltrs, frange, Us, U, p, Y, WCoefs, FrFltBnk, nh(i), nl(i)] = cwtpica(data, fs, vpo, wname, fl, fh, Peaks1, Mdl, config);
    
    H=S{1}./(6.*(std(S{1},0,2))); % normalization
    HH(i,:)=H;
    FH(i,:)=normalize(Fltrs{1},'range');
 

    L=S{2}./(15.*(std(S{2},0,2))); % normalization
    LL(i,:)=L;
    FL(i,:)=normalize(Fltrs{2},'range');
    
    audiowrite(['TempResults\CwtPicaHeart\' filesName{i}(1:end-4) '_CwtPica_Heart' '.wav'],H,fs)
    audiowrite(['TempResults\CwtPicaLung\' filesName{i}(1:end-4) '_CwtPica_Lung' '.wav'],L,fs)

end
length(find(nh))
length(find(nl))


%%
[B, ~] = sort(FH);
lqrtFH=B(ceil(.25*size(B,1)),:); 
uqrtFH=B(ceil(.75*size(B,1)),:);

[B, ~] = sort(FL);
lqrtFL=B(ceil(.25*size(B,1)),:); 
uqrtFL=B(ceil(.75*size(B,1)),:);

figure
shadedErrorBar(frange,median(FH),[uqrtFH-median(FH); -lqrtFH+median(FH)],'lineprops',{'LineWidth',1,'Color','r'});
hold on
shadedErrorBar(frange,median(FL),[uqrtFL-median(FL); -lqrtFL+median(FL)],'lineprops',{'LineWidth',1,'Color','b'});
xlabel('Hz'); ylabel('Gain');
legend('Heart filters, med \pm quart.','Lung filters, med \pm quart.','FontSize',12)
savefig('TempResults\CwtPicaFiltersShaded.fig');
saveas(gcf,'TempResults\CwtPicaFiltersShaded.eps','epsc')

figure
for i=1:size(FH,1)
    hold on
    p1=plot(frange,FH(i,:),'LineWidth',.7,'Color', [1 .75 .75 ]); 
    p2=plot(frange,FL(i,:),'LineWidth',.7,'Color',[.65 .8 1 ]);
end
hold on 
p3=plot(frange,mean(FH),'LineWidth',1.4,'Color','r');
p4=plot(frange,mean(FL),'LineWidth',1.4,'Color','b');
legend([p1 p2 p3 p4],'Heart sound filters','Lung sound filters','Mean of heart filters','Mean of lung filters','FontSize',12);
xlabel('Hz'); ylabel('Gain');
savefig('TempResults\CwtPicaFiltersOverlay.fig');
saveas(gcf,'TempResults\CwtPicaFiltersOverlay.eps','epsc')




%% Logarithmic Plots
FL(FL<=0.00001)=0.00001; FH(FH<=0.00001)=0.00001;
LogFL=20*log10(FL);
LogFH=20*log10(FH);

[B, ~] = sort(LogFH);
lqrtFH=B(ceil(.25*size(B,1)),:); 
uqrtFH=B(ceil(.75*size(B,1)),:);

[B, ~] = sort(LogFL);
lqrtFL=B(ceil(.25*size(B,1)),:); 
uqrtFL=B(ceil(.75*size(B,1)),:);

figure
shadedErrorBar(frange,median(LogFH),[uqrtFH-median(LogFH); -lqrtFH+median(LogFH)],'lineprops',{'LineWidth',1,'Color','r'});
hold on
shadedErrorBar(frange,median(LogFL),[uqrtFL-median(LogFL); -lqrtFL+median(LogFL)],'lineprops',{'LineWidth',1,'Color','b'});
xlabel('Hz'); ylabel('Gain (dB)');
legend('Heart filters, med \pm quart.','Lung filters, med \pm quart.','FontSize',12,'Location','south')
savefig('TempResults\CwtPicaLogFiltersShaded.fig');
saveas(gcf,'TempResults\CwtPicaLogFiltersShaded.eps','epsc')

figure
for i=1:size(LogFH,1)
    hold on
    p1=plot(frange,LogFH(i,:),'LineWidth',.7,'Color', [1 .75 .75 ]);
    p2=plot(frange,LogFL(i,:),'LineWidth',.7,'Color',[.65 .8 1 ]); 
end
hold on 
p3=plot(frange,mean(LogFH),'LineWidth',1.4,'Color','r');
p4=plot(frange,mean(LogFL),'LineWidth',1.4,'Color','b');
legend([p1 p2 p3 p4],'Heart sound filters','Lung sound filters','Mean of heart filters','Mean of lung filters','FontSize',12,'Location','south');
xlabel('Hz'); ylabel('Gain (dB)');
savefig('TempResults\CwtPicaLogFiltersOverlay.fig');
saveas(gcf,'TempResults\CwtPicaLogFiltersOverlay.eps','epsc')

%%
% save('TempResults\CwtPicaData.mat','FL','FH','LogFL','LogFH','HH','LL')



