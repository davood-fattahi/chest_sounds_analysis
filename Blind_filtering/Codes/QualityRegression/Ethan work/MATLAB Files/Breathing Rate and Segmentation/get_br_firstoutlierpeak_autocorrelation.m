    function heartRate=get_br_firstoutlierpeak_autocorrelation(signal_autocorrelation,max_HR,min_HR,Fs)
        min_index = round((60/max_HR)*Fs); %set max heart rate to 220
        max_index = round((60/min_HR)*Fs); %set min heart rate to 70
        
        %[peaks,locs,~,p] = findpeaks(signal_autocorrelation(min_index:max_index)); 
        [~,locs,~,p] = findpeaks(signal_autocorrelation); 
        locs=locs(locs>=min_index&locs<=max_index);
        p=p(locs>=min_index&locs<=max_index);
        if max(p)<=0.05
            bestpeak_loc=find(p==max(p),1); 
        else
            bestpeak_loc=find(p>0.05,1);
        end
        index=locs(bestpeak_loc);
        true_index = index-1;
        heartRate = 60/(true_index/Fs);
    end