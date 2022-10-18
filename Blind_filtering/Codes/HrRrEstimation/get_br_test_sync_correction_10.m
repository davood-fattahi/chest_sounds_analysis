function [BR_test]=get_br_test_sync_correction_10(Files,options,max_BR,min_BR,folder, file_ending)
BR_test=table;
for i=1:length(Files)
    file=sprintf('%ssync_%s%s.wav',folder,Files{i},file_ending);
    [audio_data,Fs]=audioread(file);
    if any(audio_data)
        BR_test(i,:)= get_br_segmentation(audio_data, Fs, max_BR,min_BR,options);
    end
end
end

