function c = correlation_evaluation(d,e)
%% Paper Information
% Hybrid Nelder-Mead search based optimal Least Mean Square algorithms for heart and lung sound separation
% https://www.sciencedirect.com/science/article/pii/S2215098616308758

%% Purpose
% To assess the performance in the extraction of desired signal

%% Inputs
% d=desired signal
% e=estimation of desired signal

% Output
% c=correlation
c=xcorr(d,e,0,'coeff');
end

