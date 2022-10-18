function r=autoCorrPeaksRatio(x,fs,AHR)



% r=autoCorrPiRatio(x, fs, AHR)
% 
% x: the signal
% fs: sampling frequency (Hz)
% AHR: average heart rate normalized by sampling frequency
%


%% Preprocessing
% Zero meaning
x=x-mean(x);

%% AutoCorrelation
Cxx=xcorr(abs(x));
Cxx=Cxx(floor(length(Cxx)./2):end);  % skipin the last samples distortion

%% envelop of AC
EnvCxx=abs(hilbert(Cxx));

%% Env peak detection based on AHR
PPeaks=find(PeakDetect(EnvCxx,AHR/fs));

r=Cxx(PPeaks(2))/Cxx(PPeaks(1));