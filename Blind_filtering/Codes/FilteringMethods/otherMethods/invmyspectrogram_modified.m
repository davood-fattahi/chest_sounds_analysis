function a = invmyspectrogram_modified(b,hop)
%% Paper Information
% Single-Channel Source Separation Tutorial Mini-Series
% https://ccrma.stanford.edu/~njb/teaching/sstutorial/ 
%% Purpose
%INVMYSPECTROGRAM Resynthesize a signal from its spectrogram.
%   See also: MYSPECTROGRAM
%% Inputs
%   B = complex array of STFT values as generated by MYSPECTROGRAM.
%   The number of rows of B is taken to be the FFT size, NFFT.
%   INVMYSPECTROGRAM resynthesizes A by inverting each frame of the 
%   FFT in B, and overlap-adding them to the output array A.  
%   NHOP is the overlap-add offset between successive IFFT frames.
%% Output
% A= time domain

[nfft,nframes] = size(b);

% nfft assumed even
No2 = nfft/2; 
a = zeros(1, nfft+(nframes-1)*hop);
 % output time offset = half of FFT size
xoff = 0 - No2;
for col = 1:nframes
  fftframe = b(:,col);
  xzp = ifft(fftframe);
  % if signal known to be real
  % xzp = real(xzp);
  x = [xzp(nfft-No2+1:nfft); xzp(1:No2)];
  % FFT's "negative-time indices" are out of range
  if xoff<0 
    ix = 1:xoff+nfft;
    % partial frames out
    a(ix) = a(ix) + x(1-xoff:nfft)'; 
  else
    ix = xoff+1:xoff+nfft;
    % overlap-add reconstruction
    a(ix) = a(ix) + x'; 
  end
  xoff = xoff + hop;
end

end