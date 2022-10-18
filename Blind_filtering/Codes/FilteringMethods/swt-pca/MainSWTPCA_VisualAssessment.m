%%%%%%%%%%%%%%%%%%%
%% Initializing
clear
close all
clc

% delete 'TempResults\*.wav'
% delete 'TempResults\*.xls'

%% Defining address of the folder containing the audio records
% DA='..\..\10sec_selected_for_separation';
% files=dir([DA '\*.mp3']); 

waitfor(InputSignals);

if strcmp(InputFormat,'all')
    files=dir([DA '\*.wav']);
elseif strcmp(InputFormat,'single')
    files.name=FileName;
end


%% setting parameters (default values)
N= 6; % number of swt levels
wname='Bior1.3'; % wavelet name
ph=[1,2]; % heart components
pl=[4,5]; % Lung components

%% configuration spread sheet
ConfigSheet=cell(size(files,1),5);
ConfigSheet(1,:)={'file name', 'number of levels in SWT', 'wavelet name', 'heart components', 'lung components'};

Save=true;
cnfrm=true;
%% do for each of the records
for i= 1 :size(files,1)
    if ~cnfrm
        break
    end    
    close all

    %%% Display the progress percentage 
    clc;disp([num2str(i) '/' num2str(size(files,1)) ', ' num2str(floor(100*i/size(files,1))) ' %']) 

    %%% Loading the signal
    [data,fs]=audioread([DA  '\' files(i).name]);    %% read the record
    
    %%% Preprocessing - Downsampling to 2000 Hz
    data=downsample(data,fs./4000); fs=4000;

    
    %%% Preprocessing - normalizing
    data=normalize(data);
  
    Repeat=true;
    while Repeat
    close all
 
%% getting parameters     
    waitfor(ParameterSetting);
    
%% SWTPCA
    [S, Fltrs, frange, Us, U, p, Y, WCoefs, FrFltBnk] =swtpca(data, fs, N, wname, (2:N) , {ph, pl} );
    
    %% plots
    psdY=normalize(pwelch(real(Y'),size(Y,2)./(.1*fs)),'range');
    figure('units','normalized','outerposition',[0.5 0 0.5 1]);
    plotcolumns(psdY);
    
    
    figure('units','normalized','outerposition',[0 0 .5 1])    
    X=real(Y);
    J=size(X,1);
    for j=1:J
        subplot(J,1,j)
        plot((1:size(X,2))./fs,normalize(X(j,:),'range'))
        xlim([0 10]); grid on; ax=gca; ax.YGrid = 'off'; set(gca,'xticklabel',{[]},'yticklabel',{[]})
        ylabel(num2str(j))
    end
    xlabel 'time (sec)'
    sgtitle 'PCA Components'
%% Component selection

    drawnow
    waitfor(ComponentSelection);   
    if ~cnfrm
        break
    end
    [S, Fltrs, frange, Us, U, p, Y, WCoefs, FrFltBnk] =swtpca(data, fs, N, wname, (2:N) , {ph, pl} );

    H=S{1}./(6.*(std(S{1},0,2))); % normalization

    L=S{2}./(15.*(std(S{2},0,2))); % normalization
%% plots

    figure(3)    
    subplot(3,1,1)
    plot((1:size(data(:),1))./fs,normalize(data(:),'range'),'Color',[.6 .2 .2])
    xlim([0 10]);    
%     set(gca,'xticklabel',(0:.5: 2.5),'FontSize',12)
    xlabel 'time (sec)'
    ylabel 'amplitutde'
    title 'raw signal'

    subplot(3,1,2)
    plot((1:length(H))./fs,normalize(H(:),'range'),'Color',[.6 .2 .2])
    xlim([0 10]);
    title 'Heart sound'
    subplot(3,1,3)
    plot((1:length(L))./fs,normalize(L(:),'range'),'Color',[.6 .2 .2])
    xlim([0 10]);
    title 'Lung sound'



    figure(4)  % frequency filtering cuased by the source estimation
    subplot(311); 

    hold on
    for k=1:size(FrFltBnk,1)
        plot(frange(1:size(FrFltBnk,2)),FrFltBnk(k,:),'LineWidth',1)
    end
    xlabel('Hz'); ylabel('amplitude'); 
    set(gca,'FontSize',12)
    title 'swt frequency subbands'

    subplot(312)
    plot(frange,normalize(Fltrs{1},'range'),'LineWidth',1); title('(b)'); xlabel('Hz'); ylabel('amplitude'); 
    set(gca,'FontSize',12)
    title 'Heart sound frequency filter'


    subplot(313)
    plot(frange,normalize(Fltrs{2},'range'),'LineWidth',1); title('(c)'); xlabel('Hz'); ylabel('amplitude'); 
    set(gca,'FontSize',12)
    title 'Lung sound frequency filter'

 %% write the results on hard disk   
    %%% Correcting the name - removing extension
    name{i}=files(i).name;
    if strcmp(name{i}(end-3:end-2),'.m') || strcmp(name{i}(end-3:end-2),'.w')
        name{i}(end-3:end)=[];
    end
    
    waitfor(Results);    
    if Save
    audiowrite(['TempResults\' name{i} '_PCAonSWT_Heart_' num2str(ph) '.wav'],H,fs)
    audiowrite(['TempResults\' name{i} '_PCAonSWT_Lung_' num2str(pl) '.wav'],L,fs)
    
%%
    ConfigSheet(i+1,:)={files(i).name, N, wname, num2str(ph), num2str(pl)};
    end
    end
end
clk=int2str(clock);  clk= clk(~isspace(clk));
writecell(ConfigSheet,['TempResults\ConfigSheet' clk '.xls']);

ConfigSheet(1,:)={'file name', 'number of levels in SWT', 'wavelet name', 'heart components', 'lung components'};

