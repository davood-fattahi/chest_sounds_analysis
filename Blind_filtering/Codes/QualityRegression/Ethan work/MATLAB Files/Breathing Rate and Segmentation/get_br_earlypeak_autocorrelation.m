    function heartRate=get_br_earlypeak_autocorrelation(signal_autocorrelation,max_HR,min_HR,Fs)
        min_index = round((60/max_HR)*Fs); %set max heart rate to 220
        max_index = round((60/min_HR)*Fs); %set min heart rate to 70
        
        %[peaks,locs,~,p] = findpeaks(signal_autocorrelation(min_index:max_index)); 
        [peaks,locs,~,p] = findpeaks(signal_autocorrelation); 
        p=p(locs>=min_index & locs<=max_index); 
        locs=locs(locs>=min_index & locs<=max_index); 
        peaks=peaks(locs>=min_index & locs<=max_index); 
        [~,I]=max(peaks);
        index=locs(I);
        val=0;
        max_p=max(p); 
        while val<0.05*Fs 
            [val,I]=min(abs(locs-index/2));
            prom=p(I); 
            if val<0.1*Fs && prom>max_p/2
                index=locs(I);
            else 
                break
            end
        end
        true_index = index-1;
        heartRate = 60/(true_index/Fs);
    end