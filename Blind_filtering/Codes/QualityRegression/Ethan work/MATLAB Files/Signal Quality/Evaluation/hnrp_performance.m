function p=hnrp_performance(HNRP)
%% Paper Information
% Blind Source Separation of Heart and Lung Sounds Based on Nonnegative Matrix Factorization
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6704646

%% Purpose
% To assess the performance in the reduction in heart noise

%% Inputs
% HNRP(%)= heart noise (or interference) reduction percentage. Which
% measures what percentage the heart sound is reduced by defining the mean
% engergy reduction
% b(%)= used percentage of HS in the mixture HLS simulation

% Output
% p(%)= performance of the evaluation

p=100*(1-abs(b-HNRP)/b);
end

