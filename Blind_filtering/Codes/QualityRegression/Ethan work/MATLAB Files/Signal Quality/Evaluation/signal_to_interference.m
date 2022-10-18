function SIR=signal_to_interference(X,M,P)
%% Paper Information
% Unsupervised SIngle-Channel Separation of Non-Stationary Signals using Gammatone Filterbank and Itakura Saito Non-negative Matrix Two-Dimensional Factorization
% https://www.researchgate.net/profile/Wai-Lok-Woo/publication/236007062_Unsupervised_Single-Channel_Separation_of_Nonstationary_Signals_Using_Gammatone_Filterbank_and_Itakura-Saito_Nonnegative_Matrix_Two-Dimensional_Factorizations/links/5efe1aed45851550508533a6/Unsupervised-Single-Channel-Separation-of-Nonstationary-Signals-Using-Gammatone-Filterbank-and-Itakura-Saito-Nonnegative-Matrix-Two-Dimensional-Factorizations.pdf

%% Purpose
% For NMF methods, to measure how well mask suppresses the interferring
% sources


%% Inputs
% X= time frequency representaiton of target
% M= time frequency mask generated in NMF method which is used in the
% mixture recording y to extract the desired target (X_opt= M.*Y)
% P= time frequency representation of interferring sources

% Output
% SIR= signal to interference ratio

SIR=norm(M.*X,'fro')^2/norm(P.*X,'fro')^2;
end