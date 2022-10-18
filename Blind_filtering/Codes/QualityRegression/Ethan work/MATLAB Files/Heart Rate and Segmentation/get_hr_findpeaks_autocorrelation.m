    function  [heartRate,locs,numpeaks]=get_hr_findpeaks_autocorrelation(signal_autocorrelation,max_HR,Fs)
        min_index = round((60/max_HR)*Fs); %set max heart rate to 220
        %max_index = round((60/min_HR)*Fs); %set min heart rate to 70
        
        [~, locs] = findpeaks(signal_autocorrelation,'MinPeakDistance' ,min_index);
        
        true_index=median(diff(locs)); 
        heartRate = 60/(true_index/Fs);
        numpeaks=length(locs); 
    end