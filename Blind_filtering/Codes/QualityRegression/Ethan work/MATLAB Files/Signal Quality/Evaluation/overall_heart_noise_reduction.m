function HNRP=overall_heart_noise_reduction(m,f)
%% Paper Information
% Blind Source Separation of Heart and Lung Sounds Based on Nonnegative Matrix Factorization
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6704646

%% Purpose
% To assess the reduction in heart noise

%% Inputs
% m= mixed heart and lung signals in time domain
% f= separated lung sound

% Output
% HNRP(%)= heart noise (or interference) reduction percentage. Which
% measures what percentage the heart sound is reduced by defining the mean
% engergy reduction

HNRP=((mean(m.^2)-mean(f.^2))/mean(m.^2))*100;
end






