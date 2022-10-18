function BR_test=get_br_test(Lung_table,BR_test,options,max_BR,min_BR)
for i=1:size(Lung_table,1)
     set=Lung_table.AnnotationSet(i);
     if set==1
         folder='/Users/ethangrooby/Downloads/Signal Quality/Annotation/Set 1 Annotation/Lung Pool/';
     elseif set==2
         folder='/Users/ethangrooby/Downloads/Signal Quality/Annotation/Set 2 Annotation/Lung Pool/';
     elseif set==3
         folder='/Users/ethangrooby/Downloads/Signal Quality/Annotation/Set 3 Annotation/Lung Pool/';
     end
     file=sprintf('%d.wav',Lung_table.index(i)); 
     [audio_data,Fs]=audioread([folder file]);
     BR_test(i,:)= get_br_segmentation(audio_data, Fs, max_BR,min_BR,options);
     %HR_value(i,1)=HR_peaks(i)*Fs*60/length(audio_data);
end
end