function tau=EstimateTimeLags(x,fs,varargin)

% tau=EstimateTimeLags(x, fs, tau0, teff, MaxLag, AHR, fh)
% 
% x: the signal
% fs: sampling frequency (Hz)
% MaxLag: maximum value for time lags (sec)
% tau0: time lags spacing (sec)
% teff: effective radius around autocorrelaton peaks. it can be scalar for 
%   effective width [tp-teff tp+teff], or vector of a pair [teffm teffp] for 
%   effective width [tp-teffm tp+teffp].
% fh: maximum frequency used for autocorrelation calculation (Hz)
% AHR: average heart rate normalized by sampling frequency
%

%% Defualt values definition
tau0=1./fs; teff= 10*tau0;  fh=(fs/2); MaxLag=2; AHR=[];


%% Getting Inputs
if nargin>=3
    tau0=varargin{1};
end
if nargin>=4
    teff=varargin{2};
end
if nargin >=5
    MaxLag=varargin{3};
end
if nargin>=6
    AHR=varargin{4};
end
if nargin==7
    fh=varargin{5};
end

if length(teff)==1
    teffm=teff; teffp=teff;
elseif length(teff)==2
    teffm=teff(1); teffp=teff(2);    
elseif length(teff)>2
    error 'wrong dimension of teff'
end

%% Preprocessing - 

% Downsampling to 2*fh
Hfr=[50 fh];
[B_H,A_H] = butter(5,[2*Hfr(1)/fs 2*Hfr(2)/fs],'bandpass');

if fh < (fs/2)
    fds=2*fh;
    x=resample(filtfilt(B_H,A_H,x),fds,fs);
    AHR=(AHR/fds);
else
    fds=fs;
    AHR=(AHR/fs);
end

% Zero meaning
x=x-mean(x);

% invering time to samples
AcMaxLag=floor(MaxLag*fds);


%% AutoCorrelation
Cxx=xcorr(x,AcMaxLag+floor(.2*fds)); % skiping the last samples distortion
Cxx=Cxx(floor(length(Cxx)./2):end);

%% envelop of AC
EnvCxx=abs(hilbert(Cxx));
EnvCxx=EnvCxx(1:AcMaxLag); % skipin the last samples distortion

%% Average Heart Rate Estimation


% AHR estimation
if isempty(AHR)
    u=uzcross(EnvCxx-mean(EnvCxx)); % upward zero crossing detection
    d=uzcross(-EnvCxx+mean(EnvCxx)); % downward zero crossing detection

    % length correcting
    u=[1;u];
    if length(u)>length(d)
        u=u(1:length(d));
    elseif length(d)>length(u)
        d=d(1:length(u));
    end
    PV=zeros(size(u));
    PI=zeros(size(u));
    for i=1:length(u)
        [PV(i), I]=max(EnvCxx(u(i):d(i)));
        PI(i)=u(i)+I-1;
    end
    [~,I]=sort(PV,'descend');
    AHR=1./(PI(I(2))-PI(I(1)));
end

%% Env peak detection based on AHR
PPeaks=floor((find(PeakDetection(EnvCxx,AHR))-1).*(fs./fds))+1;

%% proper segments finding
s=false(floor(MaxLag*fs),1);
a=PPeaks-floor(fs*teffm); b= PPeaks+floor(fs*teffp);
a(a<=0)=1; b(b>length(s))=length(s);
for i=1:length(PPeaks)
    s(a(i):b(i))=true;
end
%% time-lags finding
T=false(floor(MaxLag*fs),1);
T(1:floor(fs*tau0):floor(MaxLag*fs))=true;
tau=find(T&s);

% 
% tstmp=(1:length(Cxx))/fds;
% figure; plot(tstmp,Cxx,'LineWidth', 0.5); hold on; plot(tstmp(1:length(EnvCxx)),EnvCxx,'LineWidth', 1); 
% xlabel 'time (s)'
% ylabel amplitude
% title 'PCG autocorrelation and its envelope'
% hold on; plot(tstmp,30000*PeakDetection(EnvCxx,AHR));
% figure; hold on; plot(50.*s); plot(tau,70*ones(size(tau)),'*r');

