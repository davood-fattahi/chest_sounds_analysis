%%%%%%%%%%%%%%%%%%%
%% Initializing
clear
close all
clc
delete '*.wav'


fs=4000;
f=fs/2; % frequancy range


%% Defining address of the folder containing the audio records
DA='..\..\10sec_selected_for_separation\';
files=dir([DA '*.mp3']); 



figure('units','normalized','outerposition',[0.25 0 .4 1])    
figure('units','normalized','outerposition',[0.25 0 .4 1])    
figure('units','normalized','outerposition',[0.25 0 .4 1])    
figure('units','normalized','outerposition',[0.25 0 .4 1])    

%% do for each of the records
for i= 1 %:size(files,1)
    %%% Display the progress percentage 
    clc;disp([num2str(i) '/' num2str(size(files,1)) ', ' num2str(floor(100*i/size(files,1))) ' %']) 
    
    %%% Loading the signal
    [data,fs]=audioread([DA files(i).name]);    %% read the record
    
    %%% Preprocessing - Downsampling
    data=downsample(data,4); fs=fs/4;
    
    %%% Preprocessing - normalizing
    data=normalize(data);
    
    %%% Data properties
    l=length(data);
    
    %%% Correcting the name - removing extension
    name{i}=files(i).name;
    if name{i}(end-3:end-2)=='.m' 
        name{i}(end-3:end)=[];
    end
%% cwt
    %%% parameter setting 
    fl=50; % lower frequency limit
    fh=1000; % higher frequency limit
    vpo=4; % voice per octave
    
    %%% filterbank generation
    fb = cwtfilterbank('SignalLength',l,'Wavelet','amor', ...
    'VoicesPerOctave',vpo, 'SamplingFrequency',fs, ...
    'FrequencyLimits',[fl fh]); 
    fltrs=freqz(fb); % cwt frequency filters
    fltrs(:,1)=[];

    
    %%% cwt coefficients (complex value) 
    [Cx, f]=cwt(data,'FilterBank',fb);

    
%% pica

    %%% s1 peak detection
    Peaks=PeakDetection(abs(Cx(end-2,:)),4/fs);
    a=find(Peaks); Peaks(a(2:2:end))=0;
    figure; plot((1:size(Cx,2))./fs,abs(Cx(15,:)),'LineWidth',2); hold on; plot((1:size(Cx,2))./fs,5*Peaks); plot((1:size(Cx,2))./fs,real(Cx(end-2,:)));
    
    %%% setting the other parameters
    nw=1000; nvrlp=990; N=[]; p={(1:3),(8:14)};
    
    %%% pica
    [U,S, Cxh]=ppica_v8(Cx,Peaks,nw,nvrlp,N,p,'largestabs');
    H=icwt(Cxh{1},'amor',f,[f(end) f(1)]);
    L=icwt(Cxh{2},'amor',f,[f(end) f(1)]);
    
    
 
%% plot the sample case results 
    if i==1

        
    %%%% plot decomposed signals
    figure(1)
    W=real(Cx);
    J=size(W,1);
    for j=1:J
        subplot(J,1,j)
        plot((1:size(W,2))./fs,normalize(W(j,:),'range'))
%         xlim([5 7.5]); 
        grid on; ax=gca; ax.YGrid = 'off'; set(gca,'xticklabel',{[]},'yticklabel',{[]})
    end
%     set(gca,'xticklabel',(0:.5: 2.5),'FontSize',12)
    xlabel 'time (sec)'
    saveas(gcf,'results\cwt_dec.eps','epsc')
    
        
        
    figure(2)
    SS=real([S{1};S{2}]);
    J=size(SS,1);
    for j=1:J
        subplot(J,1,j)
        plot((1:size(SS,2))./fs,normalize(SS(j,:),'range'))
%         xlim([5 7.5]);
        grid on; ax=gca; ax.YGrid = 'off'; set(gca,'xticklabel',{[]},'yticklabel',{[]})
    end
    xlabel 'time (sec)'
%     set(gca,'xticklabel',(0:.5: 2.5),'FontSize',12)
    saveas(gcf,'results\PiCA_CWT_dec.eps','epsc')
   

%     figure(3)
%     Cxh=real(Cxh{2});
%     J=size(Cxh,1);
%     for j=1:J
%         subplot(J,1,j)
%         plot((1:size(Cxh,2))./fs,real(Cxh(j,:)))
% %         xlim([5 7.5]);
%         grid on; ax=gca; ax.YGrid = 'off'; set(gca,'xticklabel',{[]},'yticklabel',{[]})
%     end
%     xlabel 'time (sec)'
% %     set(gca,'xticklabel',(0:.5: 2.5),'FontSize',12)
%     saveas(gcf,'results\PiCA_CWT_dec.eps','epsc')
%    
       
    figure(4)    
    subplot(7,1,1)
    plot((1:size(data(:),1))./fs,normalize(data(:),'range'),'Color',[.6 .2 .2])
    xlim([5 7.5]);    set(gca,'xticklabel',(0:.5: 2.5),'FontSize',12)
    xlabel 'time (sec)'
    ylabel 'amplitutde'
    subplot(7,1,2)
    plot((1:length(H))./fs,normalize(H(:),'range'),'Color',[.6 .2 .2])
    subplot(7,1,3)
    plot((1:length(L))./fs,normalize(L(:),'range'),'Color',[.6 .2 .2])
    saveas(gcf,'results\L_H_raw_sig.eps','epsc')

        
    figure(3)       
    subplot(311); 
    ff=fs/(2*length(data)):fs/(2*length(data)):fs/2;
    hold on
    for k=1:size(fltrs,1)
        plot(ff,fltrs(k,:),'LineWidth',1)
    end
    title('(a)'); xlabel('Hz'); ylabel('amplitude'); 
    set(gca,'FontSize',10)

    subplot(312)
    im=(squeeze(sum(abs(U{1}),2)));
    image(im,'CDataMapping','scaled');
    ax=gca; ax.YTick=(1:2:length(f)); ax.YTickLabel=ceil(f(ax.YTick)); 
    tt=(1:l)./fs; ax.XTickLabel=ceil(tt(ax.XTick)); colorbar;
    title('(b)'); ylabel('frequency (Hz)'); xlabel('time (sec)'); 
    set(gca,'FontSize',10)

    
    subplot(313)
    im=(squeeze(sum(abs(U{2}),2))); im(:,[1:200 l-400:l])=0; % avoiding edge effect
    image(im,'CDataMapping','scaled');
    ax=gca; ax.YTick=(1:2:length(f)); ax.YTickLabel=ceil(f(ax.YTick)); 
    tt=(1:l)./fs; ax.XTickLabel=ceil(tt(ax.XTick)); colorbar;
    title('(c)'); ylabel('frequency (Hz)'); xlabel('time (sec)'); 
    set(gca,'FontSize',10)
    
    
    
    for k=1:size(SS,1)
        audiowrite([num2str(k) '_CWT_PiCA_sources.wav'],normalize(real(SS(k,:)))./6,fs)
    end
%     
    end



%%%% plot the sources psd    
%     window=ceil(fs/2); noverlap=ceil(fs./4); nfft=floor(size(SS,2)/5); 
%     [pxx, f] = pwelch(SS',window,noverlap,nfft,fs);
%     figure
%    for j=1:size(pxx,2)
%         subplot(size(pxx,2),1,j)
%         plot(f,pxx(:,j))
%    end


 %% write the results on hard disk   
    audiowrite([name{i} '_CWT_PiCA_Heart.wav'],normalize(H)./6,fs)
    audiowrite([name{i} '_CWT_PiCA_Lung.wav'],normalize(L)./15,fs)



end
% saveas(figure(1),'results\CWT_decomps.eps','epsc')
% saveas(figure(2),'results\after_PiCA.eps','epsc')
% saveas(figure(3),'results\raw_signal.eps','epsc')
% saveas(figure(4),'results\sample_freqfilt.eps','epsc')
% saveas(figure(5),'results\VarFreqFilts_log2.eps','epsc')
% 