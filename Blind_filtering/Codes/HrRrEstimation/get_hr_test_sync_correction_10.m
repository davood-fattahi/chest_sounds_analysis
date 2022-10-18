function [HR_test]=get_hr_test_sync_correction_10(Files,options,max_HR,min_HR,folder,file_ending)
HR_test=table;
for i=1:length(Files)
    file=sprintf('%ssync_%s%s.wav',folder,Files{i},file_ending);
    [audio_data,Fs]=audioread(file);
    if any(audio_data)
        HR_test(i,:)= get_hr_segmentation(audio_data, Fs, max_HR,min_HR,options);
    end
end
end

