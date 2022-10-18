    %% Preprocessing
    function [audio_2000_f,Fs]=get_br_preprocessing(audio_data,Fs,Fs_new)
        % dealing with padded zeros and instances of zeros
        audio_data(audio_data==0)=min(abs(audio_data(audio_data~=0))); 

        audio_2000 = resample(audio_data,Fs_new,Fs);
        Fs=Fs_new;
        
        % >150Hz 2nd order Butterworth band pass
        audio_2000_f = butterworth_high_pass_filter(audio_2000,2,150,Fs);
    end