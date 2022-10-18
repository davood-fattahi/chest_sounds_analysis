function PSR=preserved(X,M)
%% Paper Information
% Unsupervised SIngle-Channel Separation of Non-Stationary Signals using Gammatone Filterbank and Itakura Saito Non-negative Matrix Two-Dimensional Factorization
% https://www.researchgate.net/profile/Wai-Lok-Woo/publication/236007062_Unsupervised_Single-Channel_Separation_of_Nonstationary_Signals_Using_Gammatone_Filterbank_and_Itakura-Saito_Nonnegative_Matrix_Two-Dimensional_Factorizations/links/5efe1aed45851550508533a6/Unsupervised-Single-Channel-Separation-of-Nonstationary-Signals-Using-Gammatone-Filterbank-and-Itakura-Saito-Nonnegative-Matrix-Two-Dimensional-Factorizations.pdf

%% Purpose
% For NMF methods, to measure how well the generated mask for NMF method
% preserves the signal of interest x.


%% Inputs
% X= time frequency representaiton of target
% M= time frequency mask generated in NMF method which is used in the
% mixture recording y to extract the desired target (X_opt= M.*Y)

% Output
% PSR= preserved signal quality

PSR=norm(M.*X,'fro')^2/norm(X,'fro')^2;
end