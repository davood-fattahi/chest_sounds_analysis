function LSD = log_spectral_distance(d,e)
%% Paper Information
% Modulation Filtering for Heart and Lung Sound Separation from Breath Sound Recordings
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=4649547

%% Purpose
% To assess the performance in the extraction of desired signal

%% Inputs
% d= desired signal
% e= estimation of desired signal

%% Output
% LSD= average log-spectral distance between breath sound spectra and separated signal spectra

Pd= abs(fft(d)).^2/length(d);
Pe= abs(fft(e)).^2/length(e);

LSD=sqrt(mean(10*log10(Pd./Pe).^2));

end
