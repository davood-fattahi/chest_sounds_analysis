function [snr_in,snr_out] = signal_to_noise(m,e,d)
%% Paper Information
% Hybrid Nelder-Mead search based optimal Least Mean Square algorithms for heart and lung sound separation
% https://www.sciencedirect.com/science/article/pii/S2215098616308758

%% Purpose
% To assess the performance in the extraction of desired signal

%% Inputs
% d= desired signal
% e= estimation of desired signal
% m= mixture signal

%% Output
% snr_in= initial signal to noise ratio for mixture
% snr_out= signal to noise ratio after processing
snr_in=snr(m,(m-d)); %for mixture signal
snr_out=snr(e,(e-d)); %for clean signal
end



