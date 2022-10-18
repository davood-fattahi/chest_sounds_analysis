% [bSQI] = get_bSQI(audio_data,Fs,pi_vector, b_matrix,figures)
%
% A function to extract the bSQI metric, based on the agreement between two
% peak detectors within a tolerance.
%
% Implemented in:
% D. Springer et al., "Automated signal quality assessment of mobile
% phone-recorded heart sound signals," JMET, In preparation, 2016
%
%% INPUTS:
% audio_data: the PCG data
% Fs: the sampling frequency of the audio data
% pi_vector, b_matrix - trained parameters for the HMM-based method. Loaded
% from the "hmm.mat" file which was trained on different, substantial,
% database as outlined in the paper by Springer et al.
% figures: optional boolean variable dictating the display of figures
%
%% OUTPUTS:
% bSQI: the evel of agreement between two beat detectors
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

function [bSQI] = get_bSQI_modified(peak_pos_main,peak_pos_alt,Fs)
%% Find agreement between detectors:
if(isempty(peak_pos_main) || isempty(peak_pos_alt))
    F1_score  = 0;
else
    % acceptint: acceptance interval (left and right) in samples
    % As 150ms is longer than the expected duration of the S1 and S2 heart sounds, this tolerance was decreased to 100ms  
    % S1 and S2 of adults is 122 +/- 32ms and 92+/-28ms with +/- being 95%
    % CI
    % factor is 5 and 95 percentile added together for term and 17+ age group
    % Term: 120-185, 17+: 60-115
    %mean_factor=(120+185)/(60+115); 
    mean_factor=1;
    [F1_score] = Bxb_compare(peak_pos_main, peak_pos_alt, (0.1/mean_factor*Fs));
end

bSQI = F1_score;