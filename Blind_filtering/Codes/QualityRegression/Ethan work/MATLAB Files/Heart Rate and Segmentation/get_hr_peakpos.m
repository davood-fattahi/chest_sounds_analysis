   
    function [peak_pos_all,peak_pos_s1,peak_pos_s2]=get_hr_peakpos(qt,signal_length,Fs_seg,Fs)
        %% Find mid-position of qt segments
        %Split into s1 sounds, s2 sounds and FHSounds
        %S1 sounds
        s1_segments = (qt ==1);
        s2_segments = (qt ==3);

        end_points_s1 = find(diff(s1_segments));
        if(mod(length(end_points_s1),2))
            end_points_s1 = [end_points_s1, signal_length];
        end


        end_points_s2 = find(diff(s2_segments));
        if(mod(length(end_points_s2),2))

            end_points_s2 = [end_points_s2, signal_length];
        end

        mid_points_s2 = zeros(1,length(end_points_s2)/2);
        for i =1:2:length(end_points_s2)
            mid_points_s2((i+1)/2) = round((end_points_s2(i)+end_points_s2(i+1))/2);
        end


        mid_points_s1 = zeros(1,length(end_points_s1)/2);
        for i =1:2:length(end_points_s1)
            mid_points_s1((i+1)/2) = round((end_points_s1(i)+end_points_s1(i+1))/2);
        end

        %% Convert the peak positions from samples at the lower audio_segmentation_Fs to the Fs:
        peak_pos_all = round(((sort([mid_points_s1 mid_points_s2]))./Fs_seg).*Fs);
        peak_pos_s1 = round(((sort(mid_points_s1))./Fs_seg).*Fs);
        peak_pos_s2 = round(((sort(mid_points_s2))./Fs_seg).*Fs);
    end

    