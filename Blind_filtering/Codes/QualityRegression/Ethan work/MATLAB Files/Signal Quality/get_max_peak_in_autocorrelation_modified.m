% function [mSQI] = get_peak_in_autocorrelation(untruncated_autocorrelation, Fs, figures)
%
% A function to extract amplitude of the maximum peak in the autocorrelation
% signal between lag values of 0.5 and 2.33 seconds.
% Developed as part of:
% D. Springer et al., "Automated signal quality assessment of mobile
% phone-recorded heart sound signals," JMET, In preparation, 2016
% 
% This function requires the Matlab implementation of the "sampen.m" function from physionet:
% https://physionet.org/physiotools/sampen/
%
%% INPUTS:
% untruncated_autocorrelation: The untruncated autocorrelation of the audio
% data from the PCG recording. Derived in "get_autocorrelation.m".
% Fs: the sampling frequency of the autocorrelation
% figures: optional boolean variable to display figures
%
%% OUTPUTS:
% mSQI: maximum peak in autocorrelation signal between specified heart rate
% ranges
% corr_prominance: prominence value of the maximum peak based on:
% Automatic signal quality index determination of radar-recorded heart sound signals using ensemble classification
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8731709  
%% Copyright (C) 2016  David Springer
% dave.springer@gmail.com
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [mSQI, cycle_duration,corr_prominance,max_prominance,max_prominance_peak] = get_max_peak_in_autocorrelation_modified(untruncated_autocorrelation,max_HR,min_HR,Fs)


% Set the heart rate limits
% 0.43*Fs corresponds to 140 BPM
% 2*Fs corresponds to 30 BPM
%max_HR=220;
%min_HR=70; 
start_point = round((60/max_HR)*Fs);
end_point = round((60/min_HR)*Fs);

% Find the amplitude and position of the maximum peak between these two limits:
[mSQI, pos] = max(untruncated_autocorrelation(start_point:end_point));
pos1=pos;
cycle_duration = pos + start_point - 1;

%correlation prominance
% prominance value of the max peak found in mSQI
%[~,locs,~,p] = findpeaks(untruncated_autocorrelation(start_point:end_point));
%corr_prominance=p(locs==pos1);
[peaks,locs,~,p] = findpeaks(untruncated_autocorrelation);
corr_prominance=p(locs==cycle_duration);
if isempty(corr_prominance)
    corr_prominance=0;
end   
peaks=peaks(locs>=start_point & locs<=end_point); 
p=p(locs>=start_point & locs<=end_point); 
[max_prominance,I]=max(p); 
max_prominance_peak=peaks(I);  
if isempty(max_prominance_peak)
    max_prominance_peak=0;
    max_prominance=0;
end   
% if(figures)
%     figure('Name','mSQI plot');
%     plot(untruncated_autocorrelation);
%     hold on;
%     plot(pos,mSQI,'ko');
%     legend('Autocorrelation', 'Position of maximum peak between heart rate limits');
% end
end
