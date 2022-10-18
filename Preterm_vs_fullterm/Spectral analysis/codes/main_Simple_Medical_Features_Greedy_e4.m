%%%%%%%%%%%%%%%%%%%
% this function calculate the power ratio of the recordings over different
% frequency ranges, and also boxplot them and attache t-test and MWW values
% as well.


% Davood Fattahi, Feb, 2020


%% Initializing
clear
close all
clc

%%
DA='..\10sec\';
files=dir([DA '*.mp3']); %% Defining address of the folder containing the audio records


%% Extracting Features for each of the records ... 
fcl=100;  %% Low cutoff frequency
fch=2000;    %% High cutoff frequency
% w=[133 267; 200 333; 267 400; 333 467; 400 533; 467 600; 533 667; 600 733; 667 800; 733 867; 800 933; 867 999];
% w=[100 125; 125 150; 150 175; 175 200; 200 225; 225 250];
% w=[250 300; 275 325; 300 350; 325 375; 350 400; 375 425; 400 450; 425 475; 450 500];
% w=[500 600; 550 650; 600 700; 650 750; 700 800; 750 850; 800 900; 850 950; 900 1000];
% w=[400 500; 450 550; 500 600; 550 650; 600 700; 650 750; 700 800]
w=[1000 1200; 1100 1300; 1200 1400; 1300 1500; 1400 1600; 1500 1700; 1600 1800; 1700 1900; 1800 2000];
% w=[500 1000];
% w=[125 150; 450 500; 600 700];
% w=[133 267; 200 333; 267 400; 333 467];

%%%% Pre-allocation ...
Fnum=zeros(size(w,1),size(files,1));
Fname=cell(size(files,1),1);


%%%% loop for each record
for i=1:size(files,1)
    %%% Display the progress percentage 
    clc;disp([num2str(i) '/' num2str(size(files,1)) ', ' num2str(floor(100*i/size(files,1))) ' %']) 
    [data,fs]=audioread([DA files(i).name]);    %% read the record
    Fname{i}=files(i).name;
    
%    Y=nonclip(data,fs,.5,0.1,.36,'longest');
%     Y=nonclip(data,fs,1,0.1,.20,'concat');
    %%%% parameter setting for Welch
    wind = ones(1,floor(0.060*fs)); % 60 ms
    nover = floor(length(wind)*.25); % 0.25 overlap
    nfft = 2^(nextpow2(length(wind))-1); % nfft

    %%% Feature extracting
    [pxx,f] = pwelch(data,wind,nover,nfft,fs); % Welch
    for j=1:size(w,1)
            w1=w(j,1); w2=w(j,2);
            Fnum(j,i)=bandpower(pxx,f,[w1 w2],'psd')/bandpower(pxx,f,[fcl fch],'psd'); % Normalized Features
%             Fnum(j,i)=bandpower(pxx,f,[w1 w2],'psd'); % Un-Normalized Features
    end
end

%% title (name) of features
ttl=cell(size(w,1)-1,1);
ttlLog=cell(size(w,1)-1,1);
Bndttl=cell(size(w,1)-1,1);

for j=1:size(w,1)
        w1=w(j,1); w2=w(j,2);
        Bndttl(j)= cellstr([num2str(w1) '-' num2str(w2) 'Hz']);
        ttl(j)= cellstr(['P' num2str(w1) '-' num2str(w2) '/p' num2str(fcl) '-' num2str(fch) ]);
        ttlLog(j)= cellstr(['logP' num2str(w1) '-' num2str(w2) '/logp' num2str(fcl) '-' num2str(fch) ]);
end


%% loading labels

[Lnum1,Ltxt1,Lraw1]=xlsread('..\Demographics1.xlsx','Sheet1'); %% Pre terms
[Lnum2,Ltxt2,Lraw2]=xlsread('..\Demographics1.xlsx','Sheet2');  %% full terms

Lname1=Ltxt1(2:end,1);  %% Pre terms
Lname2=Ltxt2(2:end,1);  %% full terms
Lnum=[Lnum1;Lnum2];   %% all terms
Lname=[Lname1;Lname2];  %% all terms

%% removing extensions form file names
for i=1:size(Fname,1)
    if Fname{i}(end-3:end-2)=='.m' 
        Fname{i}(end-3:end)=[];
    end
end

%% Reordering the features according the labels

%finding common cases in labels and features
[ind,loc]=ismember(Lname1,Fname);

% Removing the labels not existing in features
L1=Lnum(loc~=0,:);

% Removing features not existing in labels
Fnum1=Fnum(:,loc(loc~=0));
Fnumlog1=10*log10(Fnum1);
Fname1=Fname(loc(loc~=0),:);


%finding common cases in labels and features
[ind,loc]=ismember(Lname2,Fname);

% Removing the labels not existing in features
L2=Lnum(loc~=0,:);

% Removing features not existing in labels
Fnum2=Fnum(:,loc(loc~=0));
Fnumlog2=10*log10(Fnum2);
Fname2=Fname(loc(loc~=0),:);

Fname=[Fname1; Fname2];
Fnum=[Fnum1 Fnum2];
Fnumlog=[Fnumlog1 Fnumlog2];
L=[L1; L2];
%% saving all reshaped and reordered features
Osheet=cell(size(files,1)+1,size(ttl,1)+1);
Osheet(1,1:size(ttl,1)+1)=['name' ttl'];
Osheet(2:size(files,1)+1,1)=cellstr(Fname);
Osheet(2:end,2:end)=num2cell(Fnum');
xlswrite(['..\results\PowerBandFeatures' num2str(w(1)) '-' num2str(w(end)) 'Over' num2str(fcl) '-' num2str(fch)  '.xls'],Osheet);

Osheet(1,1:size(ttlLog,1)+1)=['name' ttlLog'];
Osheet(2:size(files,1)+1,1)=cellstr(Fname);
Osheet(2:end,2:end)=num2cell(Fnumlog');
xlswrite(['..\results\LogPowerBandFeatures' num2str(w(1)) '-' num2str(w(end)) 'Over' num2str(fcl) '-' num2str(fch)  '.xls'],Osheet);


%% Correlation, Fisher measure (t-test), and MWW test calculation
%%% Pre allocation ....
CrPr=zeros(size(L,2),size(Fnum,1));
Fisher=zeros(size(Fnum,1),1);
Fisherlog=zeros(size(Fnum,1),1);
MWW=zeros(size(Fnum,1),1);
Ttest=zeros(size(Fnum,1),1);
MWWlog=zeros(size(Fnum,1),1);
Ttestlog=zeros(size(Fnum,1),1);
%%%
for i=1:size(L,2)
    for j=1:size(Fnum,1)
        CrPr(i,j)=corr(L(:,i),Fnum(j,:)','type','Pearson'); 
    end
end
for j=1:size(Fnum,1)
        Fisher(j)=((mean(squeeze(Fnum1(j,:)))-mean(squeeze(Fnum2(j,:))))^2)./(var(squeeze(Fnum1(j,:)))+var(squeeze(Fnum2(j,:))));
        Fisherlog(j )=((mean(squeeze(Fnumlog1(j,:)))-mean(squeeze(Fnumlog2(j,:))))^2)./(var(squeeze(Fnumlog1(j,:)))+var(squeeze(Fnumlog2(j,:))));
        MWW(j )= ranksum(squeeze(Fnum1(j,:)),squeeze(Fnum2(j,:)));
        [h,Ttest(j )]=ttest2(squeeze(Fnum1(j ,:)),squeeze(Fnum2(j ,:)));
        MWWlog(j )= ranksum(squeeze(Fnumlog1(j ,:)),squeeze(Fnumlog2(j ,:)));
        [h,Ttestlog(j )]=ttest2(squeeze(Fnumlog1(j ,:)),squeeze(Fnumlog2(j ,:)));
end


%% BoxPlotting 

F=Fnum'; %%% selecting 50 Hz power bandwidth
Flog=Fnumlog';  %%% selecting 50 Hz logarithm_power bandwidth

%%%% BoxPlot of powers
F1=Fnum1'; F2=Fnum2';% separating pre term and full term


position1=1:1:size(F,2);
position2=position1+0.15;
figure('units','normalized','outerposition',[0 0 1 1])
boxplot(F1,'LabelOrientation','inline','Symbol','','Color', 'b','positions', position1,'width',0.12)
set(gca,'xtick',1:size(F,2));
xt=(1:size(F,2))-.3;
yt =max([F1;F2]);
% str1=[repmat('Fisher=',size(Fisher,1),1) num2str(Fisher(:,1))];
str2=[repmat('mw=',size(MWW,1),1) num2str(MWW(:,1))];
str3=[repmat('t=',size(Ttest,1),1) num2str(Ttest(:,1))];
a=get(gca);
s=(abs(a.YLim(1)-a.YLim(2))/30);
% text(xt,yt,str1,'FontSize',14);
text(xt,yt+s,str2,'FontSize',14);
text(xt,yt+s+s,str3,'FontSize',14);

hold on
boxplot(F2,'LabelOrientation','inline','Symbol','','Labels',Bndttl(:,1),'Color', 'r','positions', position2,'width',0.12)
set(gca,'ylim',[min(min([F1;F2]),[],2) max(max([F1;F2]),[],2)+3*s]);
xlabel('frequency bands'); ylabel('Normalized Power')
saveas(gcf, ['..\results\PowerFreq' num2str(w(1)) '-' num2str(w(end)) ',NrmlzdOver' num2str(fcl) '-' num2str(fch) '.png' ])

%%%% BoxPlot of log-powers
F1=Fnumlog1'; F2=Fnumlog2';% separating pre term and full term

position1=1:1:size(Flog,2);
position2=position1+0.15;
figure('units','normalized','outerposition',[0 0 1 1])
boxplot(F1,'LabelOrientation','inline','Symbol','','Color', 'b','positions', position1,'width',0.12)
set(gca,'xtick',1:size(Flog,2));
xt=(1:size(F,2))-.3;
yt =max([F1;F2]);

% str1=[repmat('Fisher=',size(Fisherlog,1),1) num2str(Fisherlog(:,1))];
str2=[repmat('mw=',size(MWWlog,1),1) num2str(MWWlog(:,1))];
str3=[repmat('t=',size(Ttestlog,1),1) num2str(Ttestlog(:,1))];
a=get(gca);
s=(abs(a.YLim(1)-a.YLim(2))/30);
% text(xt,yt,str1,'FontSize',14);
text(xt,yt+s,str2,'FontSize',14);
text(xt,yt+s+s,str3,'FontSize',14);
hold on
boxplot(F2,'LabelOrientation','inline','Symbol','','Labels',Bndttl(:,1),'Color', 'r','positions', position2,'width',0.12)
set(gca,'ylim',[min(min([F1;F2]),[],2) max(max([F1;F2]),[],2)+3*s]);
xlabel('frequency bands'); ylabel('Normalized Power (dB)')
saveas(gcf, ['..\results\LogPowerFreq' num2str(w(1)) '-' num2str(w(end)) ',NrmlzdOver' num2str(fcl) '-' num2str(fch) '.png' ])





