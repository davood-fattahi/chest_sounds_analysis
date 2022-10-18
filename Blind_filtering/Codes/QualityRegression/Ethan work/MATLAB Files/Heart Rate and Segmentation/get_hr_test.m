function HR_test=get_hr_test(Heart_table,HR_test,options,max_HR,min_HR)
for i=1:size(Heart_table,1)
     set=Heart_table.AnnotationSet(i);
     if set==1
         folder='/Users/ethangrooby/Downloads/Signal Quality/Annotation/Set 1 Annotation/Heart Pool/';
     elseif set==2
         folder='/Users/ethangrooby/Downloads/Signal Quality/Annotation/Set 2 Annotation/Heart Pool/';
     elseif set==3
         folder='/Users/ethangrooby/Downloads/Signal Quality/Annotation/Set 3 Annotation/Heart Pool/';
     end
     file=sprintf('%d.wav',Heart_table.index(i)); 
     [audio_data,Fs]=audioread([folder file]);
     HR_test(i,:)= get_hr_segmentation(audio_data, Fs, max_HR,min_HR,options);
     %HR_value(i,1)=HR_peaks(i)*Fs*60/length(audio_data);
end
end