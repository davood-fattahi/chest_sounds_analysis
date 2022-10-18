function [FiltApprox, FiltDetail]=wavelet_filters(wname,lvl,fs)

ff=dwtfilterbank('Wavelet',wname,'Level',lvl,'SignalLength',fs);
FiltDetail=abs(freqz(ff));
FiltDetail=FiltDetail(:,floor(fs/2) + 1:end);

sf=scalingfunctions(ff);
FiltApprox=abs(fft(sf',fs));
FiltApprox=FiltApprox(1:floor(fs/2),:)';
% Filt=[FiltApprox;FiltDetail];

