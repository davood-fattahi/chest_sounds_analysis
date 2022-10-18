    function heartRate=get_hr_bestpeak_autocorrelation(signal_autocorrelation,max_HR,min_HR,Fs)
        min_index = round((60/max_HR)*Fs); %set max heart rate to 220
        max_index = round((60/min_HR)*Fs); %set min heart rate to 70

        [~,locs,~,p] = findpeaks(signal_autocorrelation(min_index:max_index)); 
        [~, bestpeak_loc]=max(p);
        index=locs(bestpeak_loc);
        true_index = index+min_index-1;

        heartRate = 60/(true_index/Fs);
    end