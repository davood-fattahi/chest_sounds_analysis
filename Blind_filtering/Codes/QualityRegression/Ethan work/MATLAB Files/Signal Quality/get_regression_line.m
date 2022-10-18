function [spectrum_slope2,intercept,r_square2]=get_regression_line(pxx,f)
%% Slope of the regression line (SL) in db/octave
meanPSD=meanfreq(pxx,f); %  The mean frequency of a power spectrum
foct=log2(f/meanPSD); % For octave we need to convert the frequency like this (taking the base as mean frequency)
spower=10*log10(pxx/pxx(end));
% To fit a linear regression line:
mdl = fitlm(foct(foct>0),spower(foct>0));
I=mdl.Coefficients{1,1};% Intercept
S=mdl.Coefficients{2,1};% slope
%limOct=log2(1000/meanPSD);

intercept=I;
spectrum_slope2 =S;
r_square2=mdl.Rsquared.Adjusted;
end
