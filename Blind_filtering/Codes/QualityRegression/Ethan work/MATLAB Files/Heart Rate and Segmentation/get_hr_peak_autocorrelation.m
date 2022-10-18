 %% Get heart rates
    function heartRate=get_hr_peak_autocorrelation(signal_autocorrelation,max_HR,min_HR,Fs)
        %% Set the max and min search indices
        % This sets the search for the highest peak in the autocorrelation to be
        % between 120 (0.5*Fs) and 30 (2*Fs) BPM
        % min_index = 0.5*Fs;
        % max_index = 2*Fs;
        %max_HR=220;
        %min_HR=70; 

        min_index = round((60/max_HR)*Fs); %set max heart rate to 220
        max_index = round((60/min_HR)*Fs); %set min heart rate to 70


        [~, index] = max(signal_autocorrelation(min_index:max_index));
        true_index = index+min_index-1;

        heartRate = 60/(true_index/Fs);
    end
    