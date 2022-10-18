function S_avg=separability_mixture(S_all)
%% Paper Information
% Unsupervised SIngle-Channel Separation of Non-Stationary Signals using Gammatone Filterbank and Itakura Saito Non-negative Matrix Two-Dimensional Factorization
% https://www.researchgate.net/profile/Wai-Lok-Woo/publication/236007062_Unsupervised_Single-Channel_Separation_of_Nonstationary_Signals_Using_Gammatone_Filterbank_and_Itakura-Saito_Nonnegative_Matrix_Two-Dimensional_Factorizations/links/5efe1aed45851550508533a6/Unsupervised-Single-Channel-Separation-of-Nonstationary-Signals-Using-Gammatone-Filterbank-and-Itakura-Saito-Nonnegative-Matrix-Two-Dimensional-Factorizations.pdf

%% Purpose
% For NMF methods, to measure the separability of overall mixture to target
% sources x1 to xk where k is total number of desired sources.


%% Inputs
% S_all= k length vector of separabilty values for individual target
% sources

% Output
% S= separability. Value of greater than 0.8 considered satisfactory
S_avg= mean(S_all);
end