function [labels, noHeart, noLung]=componentClassifyTest(x, fs, Mdl, config)



%%%% parameter setting for Welch
wind = config.wind;
nover = config.nover;
nfft = config.nfft;
fcl= config.fcl;
fch= config.fch;
fbands= config.fbands;
AHR= config.AHR;
idx=config.idx;

Features=zeros(size(x,2), size(fbands,1)+3);
[pxx,f] = pwelch(real(x),wind,nover,nfft,fs); % Welch
for i=1:size(x,2)
    bp=bandpower(pxx(:,i),f,[fcl fch],'psd');
    for j=1:size(fbands,1)
        Features(i,j)=bandpower(pxx(:,i),f,[fbands(j,1) fbands(j,2)],'psd')/bp; % Normalized Features
    end
    Peaks=find(PeakDetect(abs(x(:,i)),AHR/fs)); 
    Features(i,j+1)=std(Peaks(2:end)-Peaks(1:end-1));
    Features(i,j+2)=autoCorrPeaksRatio(x(:,i),fs,AHR);
    [~, I]=max(pxx(:,i)); Features(i,j+3)=f(I);
%     Features(i,j+4) = powerbw(pxx,f);
%     Features(i,j+5) = median(pxx); 
end
labels = predict(Mdl,normalize(Features(:,idx)));

if ~any(contains(labels,'h'))
    noHeart=true;
    [~, I]=min(abs(Features(:,8)-50));
    labels{I}='h';
else
    noHeart=false;
end
if ~any(contains(labels,'l'))
    noLung=true;
    [~, I]=min(abs(Features(:,8)-400));
    labels(I)='l';
else
    noLung=false;
end
