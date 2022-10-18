function p=finds1p(pcg,fs,hr)

% finding the s1 in the pcg, according to the peak detection of the
% signal's envelope. The envelope is calculated using absolute value of 
% hilbert transform.
% 
% inputs:
% pcg: vector of pcg data
% fs: sampling frequency
% hr: estimated heart rate per sec
% 
% output: 
% p: the peaks' position
% 
% 
% 
% Davood Fattahi, fattahi.d@gmail.com, 15/03/2020
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pcg=pcg(:); 
p=PeakDetection(abs(hilbert(pcg)),hr/fs);
end
