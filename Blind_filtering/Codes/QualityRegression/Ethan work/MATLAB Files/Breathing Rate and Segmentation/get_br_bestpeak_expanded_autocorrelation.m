function heartRate=get_br_bestpeak_expanded_autocorrelation(signal_autocorrelation,max_HR,min_HR,Fs)
min_index = round((60/max_HR)*Fs); %set max heart rate to 220
max_index = round((60/min_HR)*Fs); %set min heart rate to 70

[~,locs,~,p] = findpeaks(signal_autocorrelation);
%[~,locs,~,p] = findpeaks(signal_autocorrelation(min_index:max_index));
p=p(locs>=min_index & locs<=max_index); 
locs=locs(locs>=min_index & locs<=max_index); 
[~, bestpeak_loc]=max(p);
true_index=locs(bestpeak_loc)-1;

heartRate = 60/(true_index/Fs);
end
