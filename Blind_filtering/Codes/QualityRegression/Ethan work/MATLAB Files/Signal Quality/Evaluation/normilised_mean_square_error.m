function NMSE=normilised_mean_square_error(s,s_hat)
%% Paper Information
% On the Blind Recovery of Cardiac and Respiratory Sounds
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6879427

%% Purpose
% Assess the effectiveness in obtaining the desired signal in time domain

%% Inputs
% s= original desired signal in time domain
% s_hat= recovered desired signal in time domain
% Output
% NMSE= Normalised mean squared error in dB

NMSE=10*log10(sum((s-s_hat).^2)/sum(s.^2));
end



