% function [ccSQI] = get_ccSQI(audio_data, Fs,figures)
%
% A function to extract the cos correlation signal quality index from:
% D. Springer et al., "Automated signal quality assessment of mobile
% phone-recorded heart sound signals," JMET, In preparation, 2016
%
%% INPUTS:
% untruncated_autocorrelation: The untruncated autocorrelation of the audio
% data from the PCG recording. Derived in "get_autocorrelation.m".
% Fs: the sampling frequency of the audio recording
% figures: (optional) boolean variable for displaying figures
%
%% OUTPUTS:
% ccSQI: Correlation between the autocorrelation signal and a fitted sinusoid
%
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


function [ccSQI,HR] = get_ccSQI_modified(untruncated_autocorrelation,max_HR,min_HR, Fs)
% 120 samples either side of each peak (corresponding to the expected duration of a single heart sound for a 1kHz sampling frequency)
% factor is 5 and 95 percentile added together for term and 17+ age group
% Term: 120-185, 17+: 60-115
%mean_factor=(120+185)/(60+115);
mean_factor=1;

if(length(untruncated_autocorrelation) < 5*Fs)
    ccSQI = 0;
    HR=0;
    warning('PCG signal too short - cannot comput ccSQI');
else
    
    t = 0:1:length(untruncated_autocorrelation)-1; % Time Samples
    %sums = [];
    
    % potential range of heart rates (in cycles per second rather than
    % BPM):
    %min_HR=0.6*60; %36
    %max_HR=2.33*60; %139.8
    %max_HR=220;
    %min_HR=70;
    
    frequency_range = (min_HR/60):0.01:(max_HR/60);
    sums=zeros(1,length(frequency_range));
    for i = 1:length(frequency_range)
        
        f = frequency_range(i);% Heart beat frequency
        
        % Find the samples each side of the cosine peaks. These peaks are
        % at multiples of the "f" frequency of the cosine:
        cos_peaks1 = (1:5)'.*Fs/f - round(Fs*0.12/mean_factor);
        cos_peaks2 = (1:5)'.*Fs/f + round(Fs*0.12/mean_factor);
        count=[0;cumsum(floor(round(min(cos_peaks2,length(untruncated_autocorrelation))-min(cos_peaks1,length(untruncated_autocorrelation)),2))+1)];
        cos_spread=zeros(1,count(end));
        %cos_spread = [];
        for j = 1:5
            % Find the "spread" of the samples each side of each cosine
            % peak - cos_peaks1 is 0.12*Fs samples to the left of the peak
            % cos_peaks2 is 0.12*Fs to the right of each peak, and the
            % cos_spread finds the samples in between these two points, for
            % the five cosine peaks. This is unless the spread passes the
            % length of the signal. In that case, the spread is truncated
            % to the length.
            cos_spread(count(j)+1:count(j+1))=round(min([cos_peaks1(j) length(untruncated_autocorrelation)]):min([cos_peaks2(j) length(untruncated_autocorrelation)]));
            %cos_spread = [cos_spread, round(min([cos_peaks1(j) length(untruncated_autocorrelation)]):min([cos_peaks2(j) length(untruncated_autocorrelation)]))];
        end
        % Find the squared sum of the samples of the autocorrelation
        % function within the "cos_spread" each side of the cosine peaks
        sums(i)=sum((untruncated_autocorrelation(cos_spread)).^2); 
        %sums = [sums sum((untruncated_autocorrelation(cos_spread)).^2)];
    end
    
    % Find the heart frequency which resulted in the maximum sum of the
    % autocorrelation samples within the "cos_spread" of the fitted cosine:
    [~,b] = max(sums);
    f = frequency_range(b);
    oscil = (cos(2*pi*f/Fs*t));
    oscil = ((oscil)>0).*(oscil);
    
    % Limit the correlation to 5 seconds:
    % Do not do the above, analyze full signal
    % This equation line does not appear correct anyway (prone to errors)
    % max_cos_peaks = round(5*Fs/f) + round(Fs*0.120);
    max_cos_peaks =length(untruncated_autocorrelation); 
    correlation_coefficient = corrcoef(untruncated_autocorrelation(1:max_cos_peaks), oscil(1:max_cos_peaks));
    
    ccSQI =  correlation_coefficient(2,1);
    HR=f*60;
%     if(figures)
%         figure('Name', 'ccSQI Plot');
%         t = (1:length(untruncated_autocorrelation(1:max_cos_peaks)))./Fs;
%         plot(t,untruncated_autocorrelation(1:max_cos_peaks));
%         hold on;
%         plot(t,oscil(1:max_cos_peaks),'r');
%     end
end
end