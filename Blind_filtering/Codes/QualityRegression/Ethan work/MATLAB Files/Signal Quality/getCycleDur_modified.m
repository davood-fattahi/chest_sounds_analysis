 function dur=getCycleDur_modified(x,max_HR,min_HR,fs)
% get cardiac cycle duration from auto corelation function
% input: x, auto corelation function of an envelope, single side
%           fs, sampling frequency
%% Set the max and min search indices
% This sets the search for the highest peak in the autocorrelation to be
% between 120 (0.5*Fs) and 30 (2*Fs) BPM
% min_index = 0.5*Fs;
% max_index = 2*Fs;
%max_HR=220;
%min_HR=70; 
        
start_point = round((60/max_HR)*fs); %set max heart rate to 220
end_point = round((60/min_HR)*fs); %set min heart rate to 70
%start_point=round(0.3*fs);
%end_point=round(2*fs);

% get the maximum amplitude of the auto corelation 
[~, locs]=max(x(start_point:end_point));

dur=locs+start_point;

end