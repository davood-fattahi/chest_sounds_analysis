function E=get_spectral_energy(x,fs)
% TF representation
N=12.5/1000*fs; %12.5ms
S=0.20; %75% overlap %changed to 80% overlap due to rounding purposes
Hop_samples=round(S*N); 
window=hamming(N); 
noverlap=N-Hop_samples; 
nfft=2*N; 
% Magnitude spectrogram X
[Xcomplex,~,~]=sg(x,nfft,fs,window,noverlap); 
Xabs=abs(Xcomplex);

% spectral energy
E=sum(Xabs,1);
E = resample(E, fs, fs/(Hop_samples));
end