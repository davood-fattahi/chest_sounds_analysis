function [BR_test,BR_filt_test,Comparison]=get_br_test_sync(Sync_Data,File_Info,options,max_BR,min_BR,window_length)
count=1;
BR_test=table;
BR_filt_test=table;
Comparison=table;
for i=1:length(Sync_Data.File)
    for j=1:length(File_Info.File)
        if strcmp(Sync_Data.File{i},File_Info.File{j})
            break
        end
    end
     d_starttime= datetime(Sync_Data.Time{i}, 'Format', 'HH:mm:ss.SSSS')-seconds(window_length);
     f_starttime= datetime(File_Info.FileStartTime{j}, 'Format', 'HH:mm:ss.SSSS')+seconds(File_Info.StartTime(j));
     if d_starttime >= f_starttime 
         start= seconds(d_starttime-f_starttime); 
         if start<10-window_length
            finish= start+window_length;
            file=sprintf('/Users/ethangrooby/Downloads/Signal Quality/Annotation/Set 3 Annotation/Raw Files/Heart/sync_%s.wav',Sync_Data.File{i});
            [audio_data,Fs]=audioread(file);
            audio_data=audio_data(round(start*Fs)+1:round(finish*Fs));
            BR_test(count,:)= get_br_segmentation(audio_data, Fs, max_BR,min_BR,options);
            file=sprintf('/Users/ethangrooby/Downloads/Signal Quality/Annotation/Set 3 Annotation/Raw Files/Heart/sync_%s_filt.wav',Sync_Data.File{i});
            [audio_data,Fs]=audioread(file);
            audio_data=audio_data(round(start*Fs)+1:round(finish*Fs));
            BR_filt_test(count,:)= get_br_segmentation(audio_data, Fs, max_BR,min_BR,options);
            Comparison(count,:)=[Sync_Data(i,:) File_Info(j,4:end)]; 
            count=count+1;
         end
     end
     
end
end