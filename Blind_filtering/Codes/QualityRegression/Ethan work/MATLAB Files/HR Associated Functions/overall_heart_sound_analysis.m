function heartSound = overall_heart_sound_analysis(PCG,Fs1,window_secs)
PCG(PCG==0)=eps; %dealing with padded zeros
%
% Sample entry for the 2016 PhysioNet/CinC Challenge.
%
% INPUTS:
% recordName: string specifying the record name to process
%
% OUTPUTS:
% classifyResult: integer value where
%                     1 = abnormal recording
%                    -1 = normal recording
%                     0 = unsure (too noisy)
%
PCG=resample(PCG,4000,Fs1);
Fs1=4000; 
figures=false;
%% Load the trained parameter matrices for Springer's HSMM model.
% The parameters were trained using 409 heart sounds from MIT heart
% sound database, i.e., recordings a0001-a0409.
load('Springer_B_matrix.mat', 'Springer_B_matrix');
load('Springer_pi_vector.mat', 'Springer_pi_vector');
load('Springer_total_obs_distribution.mat', 'Springer_total_obs_distribution');

%% Load data and resample data
springer_options   = default_Springer_HSMM_options;
PCG_resampled      = resample(PCG,springer_options.audio_Fs,Fs1); % resample to springer_options.audio_Fs (1000 Hz)
%% Running runSpringerSegmentationAlgorithm.m to obtain the assigned_states
[assigned_states, heartRate, systolicTimeInterval,~,heartRate_freq,qt] = runSpringerSegmentationAlgorithm(PCG_resampled, springer_options.audio_Fs, Springer_B_matrix, Springer_pi_vector, Springer_total_obs_distribution,'Normal',false); % obtain the locations for S1, systole, s2 and diastole

[poor_quality, RMSSDall, transitCounterNorm,crit2,reducedTransitCounter]=quality_estimation_1(PCG_resampled,assigned_states);

[assigned_states_pks, ~, systolicTimeInterval_pks,heartRate_pks,~,qt_pks] = runSpringerSegmentationAlgorithm(PCG_resampled, springer_options.audio_Fs, Springer_B_matrix, Springer_pi_vector, Springer_total_obs_distribution,'Peaks',false); % obtain the locations for S1, systole, s2 and diastole
%% Running extractFeaturesFromHsIntervals.m to obtain the features for normal/abnormal heart sound classificaiton
features  = extractFeaturesFromHsIntervals(assigned_states,PCG_resampled);


%% Quality
peak_pos_1=heart_peaks(qt,springer_options,springer_options.audio_Fs,PCG_resampled);
peak_pos_2=heart_peaks(qt_pks,springer_options,springer_options.audio_Fs,PCG_resampled);

[seSQI,seSQI_1,seSQI_2,seSQI_3,seSQI_4, svdSQI, hSQI, ccSQI, sdrSQI, vSQI, mSQI, kSQI, bSQI,bSQI_1...
    ,corr_prominance,skew,Hn,max_freq,max_pow,ratio_max,coeffs,min_coeffs,max_coeffs,mean_coeffs,median_coeffs...
    ,mode_coeffs,var_coeffs,HSMM_1,HSMM_2,Ed,var_sig,var_power,f0,min_f0,max_f0,mean_f0,median_f0,mode_f0,var_f0...
    ,ZRC, output_spectral_features,periodogram_pks_features,I,S,output_lpc,output_lsf]...
    =quality_estimation_2(springer_options,figures,PCG,peak_pos_1,peak_pos_2,Fs1,assigned_states);
%% Running classifyFromHsIntervals.m to obtain the final classification result for the current recording
% not sure why I am getting errors here but not the focus on the research
% yet
%classifyResult = classifyFromHsIntervals(features.overall);
%heartSound.classifyResult = classifyResult;

%% further heart rate analysis
secs=floor(size(PCG_resampled,1)/springer_options.audio_Fs);
heartRate_seconds=zeros(1,secs);
HR_2_seconds=zeros(1,secs);
HR_2_pks_seconds=zeros(1,secs);
heartRate_pks_seconds=zeros(1,secs);
heartRate_freq_seconds=zeros(1,secs);
for k=1:secs
    if k<3
        % need sufficient size to calculate heart rate, refer to heart rate
        % calculation code
        heartRate_seconds(k)=-1;
        HR_2_seconds(k)=-1;
        HR_2_pks_seconds(k)=-1;
        heartRate_pks_seconds(k)=-1;
        heartRate_freq_seconds(k)=-1;
        continue
    elseif k<window_secs
        PCG_segment=PCG_resampled(1:k*springer_options.audio_Fs);
        assigned_states_segment=assigned_states(1:k*springer_options.audio_Fs);
        assigned_states_segment_pks=assigned_states_pks(1:k*springer_options.audio_Fs);
    else
        PCG_segment=PCG_resampled(1+(k-window_secs)*springer_options.audio_Fs:k*springer_options.audio_Fs);
        assigned_states_segment=assigned_states(1+(k-window_secs)*springer_options.audio_Fs:k*springer_options.audio_Fs);
        assigned_states_segment_pks=assigned_states_pks(1+(k-window_secs)*springer_options.audio_Fs:k*springer_options.audio_Fs);
    end
    [heartRate_seconds(k), ~,heartRate_pks_seconds(k),heartRate_freq_seconds(k)] = getHeartRateSchmidt(PCG_segment, springer_options.audio_Fs,'Normal');
    HR_2_seconds(k)=length(find(abs(diff(assigned_states_segment))>0))/4*60/(length(assigned_states_segment)/springer_options.audio_Fs);
    HR_2_pks_seconds(k)=length(find(abs(diff(assigned_states_segment_pks))>0))/4*60/(length(assigned_states_segment)/springer_options.audio_Fs);
end


%% Another option for heart rate analysis
HR_2 = length(find(abs(diff(assigned_states))>0))/4*60/(length(PCG_resampled)/springer_options.audio_Fs);
RR=diff(find(abs(diff(assigned_states))>0))/springer_options.audio_Fs; 
H2_2_alt=60/mean(RR)/4; 
HR_2_pks = length(find(abs(diff(assigned_states_pks))>0))/4*60/(length(PCG_resampled)/springer_options.audio_Fs);
RR=diff(find(abs(diff(assigned_states_pks))>0))/springer_options.audio_Fs; 
HR_2_pks_alt=60/mean(RR)/4; 

[new_peakpos] = get_Liang_peaks(PCG_resampled,springer_options.audio_Fs,figures);
HR_Liang=length(new_peakpos)/2*60/(length(PCG_resampled)/springer_options.audio_Fs);

heartSound.quality.seSQI=seSQI;
heartSound.quality.seSQI_1=seSQI_1;
heartSound.quality.seSQI_2=seSQI_2;
heartSound.quality.seSQI_3=seSQI_3;
heartSound.quality.seSQI_4=seSQI_4;
heartSound.quality.svdSQI=svdSQI;
heartSound.quality.hSQI=hSQI;
heartSound.quality.ccSQI=ccSQI;
heartSound.quality.sdrSQI=sdrSQI;
heartSound.quality.vSQI=vSQI;
heartSound.quality.mSQI=mSQI;
heartSound.quality.kSQI=kSQI;
heartSound.quality.bSQI=bSQI;
heartSound.quality.bSQI_1=bSQI_1;
heartSound.quality.overall=poor_quality;
heartSound.quality.RMSSDall=RMSSDall;
heartSound.quality.transitCounterNorm=transitCounterNorm;
heartSound.quality.crit2=crit2;
heartSound.quality.reducedTransitCounter=reducedTransitCounter;
heartSound.quality.corr_prominance=corr_prominance;
heartSound.quality.skew=skew;
heartSound.quality.Hn=Hn;
heartSound.quality.max_freq=max_freq;
heartSound.quality.max_pow=max_pow;
heartSound.quality.ratio_max=ratio_max;
heartSound.quality.coeffs=coeffs;
heartSound.quality.min_coeffs=min_coeffs;
heartSound.quality.max_coeffs=max_coeffs;
heartSound.quality.mean_coeffs=mean_coeffs;
heartSound.quality.median_coeffs=median_coeffs;
heartSound.quality.mode_coeffs=mode_coeffs;
heartSound.quality.var_coeffs=var_coeffs;
heartSound.quality.HSMM_1=HSMM_1;
heartSound.quality.HSMM_2=HSMM_2;
heartSound.quality.Ed=Ed;
heartSound.quality.var_sig=var_sig;
heartSound.quality.var_power=var_power;
heartSound.quality.f0=f0;
heartSound.quality.min_f0=min_f0;
heartSound.quality.max_f0=max_f0;
heartSound.quality.mean_f0=mean_f0;
heartSound.quality.median_f0=median_f0;
heartSound.quality.mode_f0=mode_f0;
heartSound.quality.var_f0=var_f0;
heartSound.quality.ZRC=ZRC;
heartSound.quality.output_spectral_features=output_spectral_features;
heartSound.quality.periodogram_pks_features=periodogram_pks_features;
heartSound.quality.I=I;
heartSound.quality.S=S;
heartSound.quality.output_lpc=output_lpc;
heartSound.quality.output_lsf=output_lsf;




heartSound.overall.heartRate_liang=HR_Liang;
heartSound.segments = assigned_states;
heartSound.segments_pks=assigned_states_pks;
heartSound.features = features;
heartSound.overall.heartRate = heartRate;
heartSound.overall.heartRate_seg=HR_2;
heartSound.overall.heartRate_seg_pks=HR_2_pks;
heartSound.overall.heartRate_seg_alt=H2_2_alt;
heartSound.overall.heartRate_seg_pks_alt=HR_2_pks_alt;
heartSound.overall.heartRate_pks=heartRate_pks;
heartSound.overall.heartRate_freq=heartRate_freq;
heartSound.heartRate_seconds_seg=HR_2_seconds;
heartSound.heartRate_seconds_seg_pks=HR_2_pks_seconds;
heartSound.heartRate_seconds=heartRate_seconds;
heartSound.heartRate_seconds_pks=heartRate_pks_seconds;
heartSound.heartRate_seconds_freq=heartRate_freq_seconds;
heartSound.overall.systolicTimeInterval = systolicTimeInterval;
heartSound.overall.systolicTimeInterval_pks =systolicTimeInterval_pks;

    function [assigned_states, heartRate, systolicTimeInterval,heartRate_pks,heartRate_freq,qt] = runSpringerSegmentationAlgorithm(audio_data, Fs, B_matrix, pi_vector, total_observation_distribution, option,figures)
        % function assigned_states = runSpringerSegmentationAlgorithm(audio_data, Fs, B_matrix, pi_vector, total_observation_distribution, figures)
        %
        % A function to assign states to a PCG recording using a duration dependant
        % logisitic regression-based HMM, using the trained B_matrix and pi_vector
        % trained in "trainSpringerSegmentationAlgorithm.m". Developed for use in
        % the paper:
        % D. Springer et al., "Logistic Regression-HSMM-based Heart Sound
        % Segmentation," IEEE Trans. Biomed. Eng., In Press, 2015.
        %
        %% INPUTS:
        % audio_data: The audio data from the PCG recording
        % Fs: the sampling frequency of the audio recording
        % B_matrix: the observation matrix for the HMM, trained in the
        % "trainSpringerSegmentationAlgorithm.m" function
        % pi_vector: the initial state distribution, also trained in the
        % "trainSpringerSegmentationAlgorithm.m" function
        % total_observation_distribution, the observation probabilities of all the
        % data, again, trained in trainSpringerSegmentationAlgorithm.
        % figures: (optional) boolean variable for displaying figures
        %
        %% OUTPUTS:
        % assigned_states: the array of state values assigned to the original
        % audio_data (in the original sampling frequency).
        %% Preliminary
        if(nargin < 7)
            figures = false;
        end
        
        %% Get PCG Features:
        
        [PCG_Features, featuresFs] = getSpringerPCGFeatures(audio_data, Fs);
        
        %% Get PCG heart rate
        
        [heartRate, systolicTimeInterval,heartRate_pks,heartRate_freq] = getHeartRateSchmidt(audio_data, Fs,option);
        switch option
            case 'Normal'
                [~, ~, qt] = viterbiDecodePCG_Springer(PCG_Features, pi_vector, B_matrix, total_observation_distribution, heartRate, systolicTimeInterval, featuresFs);
            case 'Peaks'
                [~, ~, qt] = viterbiDecodePCG_Springer(PCG_Features, pi_vector, B_matrix, total_observation_distribution, heartRate_pks, systolicTimeInterval, featuresFs);
            otherwise
                [~, ~, qt] = viterbiDecodePCG_Springer(PCG_Features, pi_vector, B_matrix, total_observation_distribution, heartRate, systolicTimeInterval, featuresFs);
        end
        
        
        assigned_states = expand_qt(qt, featuresFs, Fs, length(audio_data));
        
        if(figures)
            figure('Name','Derived state sequence');
            t1 = (1:length(audio_data))./Fs;
            plot(t1,normalise_signal(audio_data),'k');
            hold on;
            plot(t1,assigned_states,'r--');
            legend('Audio data', 'Derived states');
        end
        
    end

    function features = extractFeaturesFromHsIntervals(assigned_states, PCG)
        %
        % This function calculate 20 features based on the assigned_states by running "runSpringerSegmentationAlgorithm.m" function
        %
        % INPUTS:
        % assigned_states: the array of state values assigned to the sound recording.
        % PCG: resampled sound recording with 1000 Hz.
        %
        % OUTPUTS:
        % features: the obtained 20 features for the current sound recording
        %
        %
        % Written by: Chengyu Liu, January 22 2016
        %             chengyu.liu@emory.edu
        %
        % Last modified by:
        %
        %
        % $$$$$$ IMPORTANT
        % Please note: the calculated 20 features are only some pilot features, some features maybe
        % helpful for classifying normal/abnormal heart sounds, some maybe
        % not. You need re-construct the features for a more accurate classification.
        
        
        %% We just assume that the assigned_states cover at least 2 whole heart beat cycle
        indx = find(abs(diff(assigned_states))>0); % find the locations with changed states
        
        if assigned_states(1)>0   % for some recordings, there are state zeros at the beginning of assigned_states
            switch assigned_states(1)
                case 4
                    K=1;
                case 3
                    K=2;
                case 2
                    K=3;
                case 1
                    K=4;
            end
        else
            switch assigned_states(indx(1)+1)
                case 4
                    K=1;
                case 3
                    K=2;
                case 2
                    K=3;
                case 1
                    K=0;
            end
            K=K+1;
        end
        
        indx2                = indx(K:end);
        rem                  = mod(length(indx2),4);
        indx2(end-rem+1:end) = [];
        A                    = reshape(indx2,4,length(indx2)/4)'; % A is N*4 matrix, the 4 columns save the beginnings of S1, systole, S2 and diastole in the same heart cycle respectively
        
        %% Feature calculation
        m_RR        = round(mean(diff(A(:,1))));             % mean value of RR intervals
        sd_RR       = round(std(diff(A(:,1))));              % standard deviation (SD) value of RR intervals
        mean_IntS1  = round(mean(A(:,2)-A(:,1)));            % mean value of S1 intervals
        sd_IntS1    = round(std(A(:,2)-A(:,1)));             % SD value of S1 intervals
        mean_IntS2  = round(mean(A(:,4)-A(:,3)));            % mean value of S2 intervals
        sd_IntS2    = round(std(A(:,4)-A(:,3)));             % SD value of S2 intervals
        mean_IntSys = round(mean(A(:,3)-A(:,2)));            % mean value of systole intervals
        sd_IntSys   = round(std(A(:,3)-A(:,2)));             % SD value of systole intervals
        mean_IntDia = round(mean(A(2:end,1)-A(1:end-1,4)));  % mean value of diastole intervals
        sd_IntDia   = round(std(A(2:end,1)-A(1:end-1,4)));   % SD value of diastole intervals
        
        col=size(A,1)-1;
        R_SysRR= zeros(1,col);
        R_DiaRR= zeros(1,col);
        R_SysDia= zeros(1,col);
        P_S1= zeros(1,col);
        P_Sys= zeros(1,col);
        P_S2= zeros(1,col);
        P_Dia= zeros(1,col);
        P_SysS1= zeros(1,col);
        P_DiaS2= zeros(1,col);
        
        for i=1:size(A,1)-1
            R_SysRR(i)  = (A(i,3)-A(i,2))/(A(i+1,1)-A(i,1))*100;
            R_DiaRR(i)  = (A(i+1,1)-A(i,4))/(A(i+1,1)-A(i,1))*100;
            R_SysDia(i) = R_SysRR(i)/R_DiaRR(i)*100;
            
            P_S1(i)     = sum(abs(PCG(A(i,1):A(i,2))))/(A(i,2)-A(i,1));
            P_Sys(i)    = sum(abs(PCG(A(i,2):A(i,3))))/(A(i,3)-A(i,2));
            P_S2(i)     = sum(abs(PCG(A(i,3):A(i,4))))/(A(i,4)-A(i,3));
            P_Dia(i)    = sum(abs(PCG(A(i,4):A(i+1,1))))/(A(i+1,1)-A(i,4));
            if P_S1(i)>0
                P_SysS1(i) = P_Sys(i)/P_S1(i)*100;
            else
                P_SysS1(i) = 0;
            end
            if P_S2(i)>0
                P_DiaS2(i) = P_Dia(i)/P_S2(i)*100;
            else
                P_DiaS2(i) = 0;
            end
        end
        
        m_Ratio_SysRR   = mean(R_SysRR);  % mean value of the interval ratios between systole and RR in each heart beat
        sd_Ratio_SysRR  = std(R_SysRR);   % SD value of the interval ratios between systole and RR in each heart beat
        m_Ratio_DiaRR   = mean(R_DiaRR);  % mean value of the interval ratios between diastole and RR in each heart beat
        sd_Ratio_DiaRR  = std(R_DiaRR);   % SD value of the interval ratios between diastole and RR in each heart beat
        m_Ratio_SysDia  = mean(R_SysDia); % mean value of the interval ratios between systole and diastole in each heart beat
        sd_Ratio_SysDia = std(R_SysDia);  % SD value of the interval ratios between systole and diastole in each heart beat
        
        indx_sys = find(P_SysS1>0 & P_SysS1<100);   % avoid the flat line signal
        if length(indx_sys)>1
            m_Amp_SysS1  = mean(P_SysS1(indx_sys)); % mean value of the mean absolute amplitude ratios between systole period and S1 period in each heart beat
            sd_Amp_SysS1 = std(P_SysS1(indx_sys));  % SD value of the mean absolute amplitude ratios between systole period and S1 period in each heart beat
        else
            m_Amp_SysS1  = 0;
            sd_Amp_SysS1 = 0;
        end
        indx_dia = find(P_DiaS2>0 & P_DiaS2<100);
        if length(indx_dia)>1
            m_Amp_DiaS2  = mean(P_DiaS2(indx_dia)); % mean value of the mean absolute amplitude ratios between diastole period and S2 period in each heart beat
            sd_Amp_DiaS2 = std(P_DiaS2(indx_dia));  % SD value of the mean absolute amplitude ratios between diastole period and S2 period in each heart beat
        else
            m_Amp_DiaS2  = 0;
            sd_Amp_DiaS2 = 0;
        end
        
        features.general.m_RR = m_RR;
        features.general.sd_RR = sd_RR;
        features.general.mean_IntS1 = mean_IntS1;
        features.general.sd_IntS1 = sd_IntS1;
        features.general.mean_IntS2 = mean_IntS2;
        features.general.sd_IntS2 = sd_IntS2;
        features.general.mean_IntSys = mean_IntSys;
        features.general.sd_IntSys = sd_IntSys;
        features.general.mean_IntDia = mean_IntDia;
        features.general.sd_IntDia = sd_IntDia;
        features.ratio.m_Ratio_SysRR = m_Ratio_SysRR;
        features.ratio.sd_Ratio_SysRR = sd_Ratio_SysRR;
        features.ratio.m_Ratio_DiaRR = m_Ratio_DiaRR;
        features.ratio.sd_Ratio_DiaRR = sd_Ratio_DiaRR;
        features.ratio.m_Ratio_SysDia = m_Ratio_SysDia;
        features.ratio.sd_Ratio_SysDia = sd_Ratio_SysDia;
        features.amp.m_Amp_SysS1 = m_Amp_SysS1;
        features.amp.sd_Amp_SysS1 = sd_Amp_SysS1;
        features.amp.m_Amp_DiaS2 = m_Amp_DiaS2;
        features.amp.sd_Amp_DiaS2 = sd_Amp_DiaS2;
        
        features.overall = [m_RR sd_RR  mean_IntS1 sd_IntS1  mean_IntS2...
            sd_IntS2  mean_IntSys sd_IntSys  mean_IntDia sd_IntDia...
            m_Ratio_SysRR sd_Ratio_SysRR m_Ratio_DiaRR sd_Ratio_DiaRR...
            m_Ratio_SysDia sd_Ratio_SysDia m_Amp_SysS1 sd_Amp_SysS1...
            m_Amp_DiaS2 sd_Amp_DiaS2];
    end



% function [heartRate systolicTimeInterval] = getHeartRateSchmidt(audio_data, Fs, figures)
%
% Derive the heart rate and the sytolic time interval from a PCG recording.
% This is used in the duration-dependant HMM-based segmentation of the PCG
% recording.
%
% This method is based on analysis of the autocorrelation function, and the
% positions of the peaks therein.
%
% This code is derived from the paper:
% S. E. Schmidt et al., "Segmentation of heart sound recordings by a
% duration-dependent hidden Markov model," Physiol. Meas., vol. 31,
% no. 4, pp. 513-29, Apr. 2010.
%
% Developed by David Springer for comparison purposes in the paper:
% D. Springer et al., "Logistic Regression-HSMM-based Heart Sound
% Segmentation," IEEE Trans. Biomed. Eng., In Press, 2015.
%
%% INPUTS:
% audio_data: The raw audio data from the PCG recording
% Fs: the sampling frequency of the audio recording
% figures: optional boolean to display figures
%
%% OUTPUTS:
% heartRate: the heart rate of the PCG in beats per minute
% systolicTimeInterval: the duration of systole, as derived from the
% autocorrelation function, in seconds
%
%% Copyright (C) 2016  David Springer
% dave.springer@gmail.com
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

    function [heartRate, systolicTimeInterval, heartRate_pks,heartRate_freq] = getHeartRateSchmidt(audio_data, Fs,option ,figures)
        
        if nargin < 4
            figures = false;
        end
        
        %% Get heatrate:
        % From Schmidt:
        % "The duration of the heart cycle is estimated as the time from lag zero
        % to the highest peaks between 500 and 2000 ms in the resulting
        % autocorrelation"
        % This is performed after filtering and spike removal:
        
        %% 25-400Hz 4th order Butterworth band pass
        audio_data = butterworth_low_pass_filter(audio_data,2,400,Fs, false);
        audio_data = butterworth_high_pass_filter(audio_data,2,25,Fs);
        
        %% Spike removal from the original paper:
        audio_data = schmidt_spike_removal(audio_data,Fs);
        
        %% Find the homomorphic envelope
        homomorphic_envelope = Homomorphic_Envelope_with_Hilbert(audio_data, Fs);
        
        %% Find the autocorrelation:
        y=homomorphic_envelope-mean(homomorphic_envelope);
        [c] = xcorr(y,'coeff');
        signal_autocorrelation = c(length(homomorphic_envelope)+1:end);
        
        %min_index = 0.5*Fs; %adult max heart arate set to 120
        %max_index = 2*Fs; %adult min heart rate set to 30
        
        
        max_HR=220;
        min_HR=70; %70
        %max frequency
        Y = fft(y);
        L=length(y);
        P2 = abs(Y/L);
        P1 = P2(1:round(L/2)+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = Fs*(0:(L/2))/L;
        reduced_f=f(f>min_HR/60 & f<max_HR/60);
        reduced_P1=P1(f>min_HR/60 & f<max_HR/60);
        heartRate_freq=reduced_f(reduced_P1==max(reduced_P1))*60;
        
        
        min_index = round((60/max_HR)*Fs); %set max heart rate to 220
        max_index = round((60/min_HR)*Fs); %set min heart rate to 70
        
        [~, index] = max(signal_autocorrelation(min_index:max_index));
        true_index = index+min_index-1;
        
        heartRate = 60/(true_index/Fs);
        [~, locs] = findpeaks(y,'MinPeakDistance' ,Fs/4,'MinPeakHeight',median(y));
        duration=length(y)/Fs;
        %heartRate_pks=length(locs)/duration*60;
        heartRate_pks=60/(mean(diff(locs))/Fs);
        %x=0:1/fs:(duration-1/fs);
        %plot(x,y,x(locs),y(locs),'r*')
        %plot(x,audio_data,x(locs),audio_data(locs),'r*')
        
        %% Find the systolic time interval:
        % From Schmidt: "The systolic duration is defined as the time from lag zero
        % to the highest peak in the interval between 200 ms and half of the heart
        % cycle duration"
        
        switch option
            case 'Normal'
                max_sys_duration = round(((60/heartRate)*Fs)/2);
            case 'Peaks'
                max_sys_duration = round(((60/heartRate_pks)*Fs)/2);
            otherwise
                max_sys_duration = round(((60/heartRate)*Fs)/2);
        end
        
        min_sys_duration = round(0.2*Fs);  %currently corresponds to 150bpm. should edit to 250 which is (60/250)*Fs)/2
        if max_sys_duration<=min_sys_duration
            max_sys_duration=min_sys_duration+1;
        end
        
        [~, pos] = max(signal_autocorrelation(min_sys_duration:max_sys_duration));
        systolicTimeInterval = (min_sys_duration+pos)/Fs;
        
        
        if(figures)
            figure('Name', 'Heart rate calculation figure');
            plot(signal_autocorrelation);
            hold on;
            plot(true_index, signal_autocorrelation(true_index),'ro');
            plot((min_sys_duration+pos), signal_autocorrelation((min_sys_duration+pos)), 'mo');
            xlabel('Samples');
            legend('Autocorrelation', 'Position of max peak used to calculate HR', 'Position of max peak within systolic interval');
        end
    end

    function [poor_quality, RMSSDall, transitCounterNorm,crit2,reducedTransitCounter]=quality_estimation_1(PCG_resampled,assigned_states)
        % PCG Classification Using a Neural Network Approach by Iga
        % Grzegorczyk et al.
        %Assessment of signal quality
        % We assessed the quality of the PCG signal as a first,
        % starting step. At the beginning we used the method for
        % preliminary detection of peaks in signal, which are
        % inclined to occur during the S1 or S2 phase (whichever
        % have higher amplitude in PCG cycle). We calculated the
        % energy of the normalized signal (the normalization was
        % performed separately in each of 25 windows of the
        % signal) as described in [4]. Next we determined the peaks
        % with Giera?towski et al. method [5] based on the slope
        % detection in the signal. In the next part of the method, we
        % calculate Wavelet coefficients of the whole signal by
        % using Daubechies-2 wavelet at second decomposition
        % level (according to [6]) and use these coefficients as an
        % input signal to later evaluations. We created three sets of
        % criteria for assessing signal quality. If the signal did not
        % meet at least one criterion it was classified as uncertain.
        % First criterion: RMSSD (root mean square of successive
        % differences) of the signal must be lower or equal to 0.026
        % [arbitrary units]. Moreover, the number of times which
        % the signal intersected the horizontal line determined by
        % 0.85 quantile of values of that signal divided by the signal
        % length must be lower than 0.06. Second criterion: we
        % analyzed the signal in 2200 ms length moving windows
        % with overlap 25% and checked how many peaks are
        % detected in each window. We assigned 1 to each window
        % containing 2-4 peaks and 0 in other cases. The criterion
        % was met when the number of scores equal to 1 is 65% of
        % the all notes calculated in each window. Third criterion
        % was similar to the first one, but in this case there was
        % calculated the value of 0.58 percentile of signal and the
        % number of detected peaks was extracted from the whole
        % number of intersections. According to these rules, the
        % number of intersections must be fewer than 18. Note that
        % the second and the third criterion were checked
        % depending on the number of detected peaks per one
        % second of the signal. If that value was higher than 1.1 s,
        % signal had to meet only the second criterion. Otherwise,
        % only the third criterion was checked.
        
        
        %Criteria 1: RMSSDall transitCounterNorm transitCounterNorm>0.06 && RMSSDall>=0.026
        %Criteria 2: crit2<=0.65
        %Criteria 3: reducedTransitCounter>=18
        [~, ~, ~, ~, PCG_Hi_Hi] = falkowyFeature(PCG_resampled,assigned_states);
        [ qrs_pos,  PCG_Hi_Hi_seg2]=checkingZero2(PCG_Hi_Hi);
        [ifzero, RMSSDall, transitCounterNorm,crit2,reducedTransitCounter]=detectNoisySignals2(PCG_Hi_Hi_seg2, qrs_pos);
        poor_quality=ifzero; %ifzero==1 -> classifyResult=0; i.e. poor quality
    end
    function [seSQI,seSQI_1,seSQI_2,seSQI_3,seSQI_4, svdSQI, hSQI, ccSQI, sdrSQI, vSQI, mSQI, kSQI, bSQI,bSQI_1...
            ,corr_prominance,skew,Hn,max_freq,max_pow,ratio_max,coeffs,min_coeffs,max_coeffs,mean_coeffs,median_coeffs...
            ,mode_coeffs,var_coeffs,HSMM_1,HSMM_2,Ed,var_sig,var_power,f0,min_f0,max_f0,mean_f0,median_f0,mode_f0,var_f0...
            ,ZRC, output_spectral_features,periodogram_pks_features,I,S,output_lpc,output_lsf]...
            =quality_estimation_2(springer_options,figures,audio_data,peaks_1,peaks_2,Fs1,assigned_states)
        %% Test feature extraction:
        % Load the options
        % Load the hidden Markov model parameters that have been trained on a
        % substantial database and saved. This loads the pi_vector and b_matrix
        % variables which are passed to the getSignalQualityIndices.m function:
        load('hmm.mat','b_matrix','pi_vector');
        % Extract the features
        [seSQI,seSQI_1,seSQI_2,seSQI_3,seSQI_4, svdSQI, hSQI, ccSQI, sdrSQI, vSQI, mSQI, kSQI, bSQI,bSQI_1...
            ,corr_prominance,skew,Hn,max_freq,max_pow,ratio_max,coeffs,min_coeffs,max_coeffs,mean_coeffs,median_coeffs...
            ,mode_coeffs,var_coeffs,HSMM_1,HSMM_2,Ed,var_sig,var_power,f0,min_f0,max_f0,mean_f0,median_f0,mode_f0,var_f0...
            ,ZRC, output_spectral_features,periodogram_pks_features,I,S,output_lpc,output_lsf]...
            = getSignalQualityIndices(audio_data, Fs1,figures,springer_options,peaks_1,peaks_2,pi_vector, b_matrix,assigned_states);
    end

    function [seSQI,seSQI_1,seSQI_2,seSQI_3,seSQI_4, svdSQI, hSQI, ccSQI, sdrSQI, vSQI, mSQI, kSQI, bSQI,bSQI_1...
            ,corr_prominance,skew,Hn,max_freq,max_pow,ratio_max,coeffs,min_coeffs,max_coeffs,mean_coeffs,median_coeffs...
            ,mode_coeffs,var_coeffs,HSMM_1,HSMM_2,Ed,var_sig,var_power,f0,min_f0,max_f0,mean_f0,median_f0,mode_f0,var_f0...
            ,ZRC, output_spectral_features,periodogram_pks_features,I,S,output_lpc,output_lsf]...
            = getSignalQualityIndices(audio_data, Fs,figures,springer_options,peaks_1,peaks_2,pi_vector, b_matrix,assigned_states)
        % function [seSQI, svdSQI, hSQI, ccSQI, sdrSQI, vSQI, mSQI, kSQI, bSQI] = getSignalQualityIndices(audio_data, Fs,pi_vector, b_matrix,figures)
        %
        % A function to extract signal quality indices from heart sound recordings.
        % These indices are based on the paper:
        % D. Springer et al., "Automated signal quality assessment of mobile
        % phone-recorded heart sound signals," JMET, In preparation, 2016
        %
        %% INPUTS:
        % audio_data: The audio data from the PCG recording
        % Fs: the sampling frequency of the audio recording
        % pi_vector, b_matrix: the HMM parameters needed for the bSQI. Instead of
        % being loaded for each signal, these should be loaded from the "hmm.mat"
        % file once and passed to this function.
        % figures: (optional) boolean variable for displaying figures
        %
        %% OUTPUTS:
        % [seSQI, svdSQI, hSQI, ccSQI, sdrSQI, vSQI, mSQI, kSQI, bSQI]
        % These are the nine signal quality indices derived from the PCG recordings
        % as outlined in the above publication
        %
        %% Copyright (C) 2016  David Springer
        % dave.springer@gmail.com
        %
        % This program is free software: you can redistribute it and/or modify
        % it under the terms of the GNU General Public License as published by
        % the Free Software Foundation, either version 3 of the License, or
        % any later version.
        %
        % This program is distributed in the hope that it will be useful,
        % but WITHOUT ANY WARRANTY; without even the implied warranty of
        % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        % GNU General Public License for more details.
        %
        % You should have received a copy of the GNU General Public License
        % along with this program.  If not, see <http://www.gnu.org/licenses/>.
        %% Prelim
        if nargin < 3
            figures = false;
        end
        %springer_options = default_Springer_Signal_Quality_options;
        
        %% Resample to the frequency set in the options file:
        audio_downsampled = resample(audio_data,springer_options.audio_Fs,Fs);
        
        %% Remove noise spikes in the signal using the method developed by Schmidt et al:
        [audio_downsampled] = schmidt_spike_removal(audio_downsampled,springer_options.audio_Fs);
        
        %% Get the truncated and untrancated autocorrelation functions:
        [truncated_autocorrelation, untruncated_autocorrelation] = get_autocorrelation(audio_downsampled,springer_options.audio_Fs,figures);
        semi_truncated_autocorrelation=untruncated_autocorrelation(find(untruncated_autocorrelation==truncated_autocorrelation(1),1):end); 
        max_HR=220;
        min_HR=70;
        %% Get the sample_entropy sqi:
        % Signal Quality Classification of mobile phone-recorded
        % phonocardiogram signals 2014 David Springer et al.
        %M=1,r=0.01
        %M=2, r=0.01
        %M=1, r=0.001
        %M=2, r=0.001
        M = 2;
        r = 0.0008;
        %seSQI = get_entropy(truncated_autocorrelation, M, r);
        seSQI = get_entropy(semi_truncated_autocorrelation, M, r);
        M=1;
        r=0.01;
        %seSQI_1= get_entropy(truncated_autocorrelation, M, r);
        seSQI_1= get_entropy(semi_truncated_autocorrelation, M, r);
        M=2;
        r=0.01;
        %seSQI_2= get_entropy(truncated_autocorrelation, M, r);
        seSQI_2= get_entropy(semi_truncated_autocorrelation, M, r);
        M=1;
        r=0.001;
        %seSQI_3= get_entropy(truncated_autocorrelation, M, r);
        seSQI_3= get_entropy(semi_truncated_autocorrelation, M, r);
        M=2;
        r=0.001;
        %seSQI_4= get_entropy(truncated_autocorrelation, M, r);
        seSQI_4= get_entropy(semi_truncated_autocorrelation, M, r);
        
        
        %% Get the SVD sqi:
        %[svdSQI] = get_SVD_score(truncated_autocorrelation,springer_options.audio_Fs,max_HR,min_HR);
        [svdSQI] = get_SVD_score(semi_truncated_autocorrelation,springer_options.audio_Fs,max_HR,min_HR);
        
        %% Get Hjorth activity score
        %[hSQI] = get_hjorth_activity_score(truncated_autocorrelation);
        [hSQI] = get_hjorth_activity_score(semi_truncated_autocorrelation);
        
        %% get ccSQI:
        [ccSQI] = get_ccSQI(untruncated_autocorrelation, springer_options.audio_Fs,max_HR,min_HR, figures);
        
        %% Get the sdrSQI
        % Changed to 4000Hz for lung sound analysis
        audio_downsampled_2000 = resample(audio_data,4000,Fs);
        [sdrSQI] = get_signal_to_noise_psd_score(audio_downsampled_2000,4000);
        
        %% Get the vSQI
        %[vSQI] = get_variance_score(truncated_autocorrelation);
        [vSQI] = get_variance_score(semi_truncated_autocorrelation);
        
        %% Get mSQI:
        [mSQI, ~,corr_prominance] = get_max_peak_in_autocorrelation(untruncated_autocorrelation, springer_options.audio_Fs,max_HR,min_HR,figures);
        
        %% get kSQI:
        [kSQI] = get_kurtosis_score(audio_data);
        
        %% get bSQI
        [bSQI] =get_bSQI(audio_data,Fs,pi_vector, b_matrix,figures);
        [bSQI_1] = get_bSQI_1(peaks_1,peaks_2,Fs);
        
        
        %Automatic Signal Quality Index Determination of Radar-Recorded Heart Sound Signals Using Ensemble Classification
        %Kilin Shi et al. 2020
        %skewness
        skew = skewness(audio_data);
        
        %spectral entropy
        %S(m)= abs square for fft
        %pm=S(m)/sum(S)
        %Hn=-(sum(pm*logpm))/log2N
        %N= total number of frequency points
        NFFT=1028;
        spec=fft(audio_data,NFFT);
        S=abs(spec(1:NFFT/2+1));
        N=length(S);
        p=S/sum(S);
        Hn=-(sum(p.*log(p)))/log2(N);
        
        
        %dominant frequency features
        %max freq, its value and ratio of max/total energy
        NFFT=1028;
        %L = length(audio_data);   % Length of signal
        f = Fs*(0:(NFFT/2))/NFFT;
        spec=fft(audio_data,NFFT);
        pow=2*abs(spec(1:NFFT/2+1));
        
        max_freq=f(pow==max(pow));
        max_pow=max(pow);
        ratio_max=max(pow)/sum(pow);
        
        
        %mfcc
        %13 features
        %Frame the window using a hamming window and a window length of 25 ms with a sliding window of 10 ms.
        coeffs = mfcc(audio_data,Fs,'WindowLength',round(0.025*Fs),...
            'OverlapLength',round(0.015*Fs),'NumCoeffs',13);
        min_coeffs=min(coeffs);
        max_coeffs=max(coeffs);
        mean_coeffs=mean(coeffs);
        median_coeffs=median(coeffs);
        mode_coeffs=mode(coeffs);
        var_coeffs=var(coeffs);
        %That is, the min, max, mean, median, mode, and variance of each coefficient in each frame
        % need to check this one but is done at each level i.e. 13 *6
        % values
        coeffs = mfcc(audio_data,Fs,'WindowLength',length(audio_data),'NumCoeffs',13);
        
        
        %correlation prominance
        % prominance value of the max peak found in mSQI
        %corr_prominance
        
        %HSMM Quality Factor (HSMM-QF)
        %heart sound segmentation
        % dividing the mean height of the envelope during S1/S2 by the
        % average mean height of diastole and systole. Two features, one
        % for each heart sound, result from this. The correctness of the
        % segmentation is not crucial since a wrong segmentation, e.g.,
        % in case of noise, also indicates a low signal quality which in
        % that case would lead to a low HSMM-QF value. Since in this case
        % the envelope height during S1/S2 would roughly be equal to that
        % during the diastole/systole, the HSMM-QF would be approximately
        % 1. High-quality heart sounds signals however are expected to
        % result in a high HSMM-QF values ? 1.
        %assigned_states
        %untruncated_autocorrelation
        locS1=assigned_states==1;
        locSys=assigned_states==2;
        locS2=assigned_states==3;
        locDia=assigned_states==4;
        locSysDia= (locSys|locDia);
        HSMM_1=abs(mean(untruncated_autocorrelation(locS1))/mean(untruncated_autocorrelation(locSysDia)));
        HSMM_2=abs(mean(untruncated_autocorrelation(locS2))/mean(untruncated_autocorrelation(locSysDia)));
        
        %Improving the Quality of Point of Care Diagnostics with Real-Time Machine Learning in Low Literacy LMIC Settings
        
        %Wavelet features Using these frequency ranges as reference, the
        %DUS recordings were de- composed into 4 levels using the Discrete
        %Wavelet multi-resolution analysis (Table 2). The reverse
        %biorthogonal wavelet rbio3.9 was selected as the mother wavelet
        %since it was able to correctly classify more good and poor quality
        %segments than other mother wavelets for DUS recordings made with
        %the same transducer used in this study [36, 38].
        %As a feature, the percentage of energy content at each decompo- sition level was computed as follows:
        wname='rbio3.9';
        N=4;
        %Discrete stationary wavelet transform 1-D
        % correcting the length of signal
        fin=size(audio_data,1)-mod(size(audio_data,1),2^N);
        data_e1=audio_data(1:fin);
        [SWA,SWD]=swt(data_e1,N,wname);
        SWC= real([SWA ; SWD]);
        Et=sum(sum(SWC.^2));
        Ed=100*sum(SWC'.^2)/Et;
        
        
        %Statistics
        % The first SQI was the variance, measuring total power or by how the signal varied about its mean.
        var_sig=var(audio_data);
        var_power=var(pow);
        
        %25ms window Estimate fundamental frequency hamming window and a sliding window
        %of 10 ms (15 ms overlap), which are the most common window?s
        %length for speech recognition [13]. For each frame, the F 0 was
        %estimated by taking the inverse of the quefrency within 50 and
        %1000 Hz with the highest peak. After the F 0 was estimated for
        %each possible window of the segment, all the values were
        %distributed in a histogram, calculating the number of bins from
        %the average provided by Sturgis? [40] and Rice?s [16] rules. The F
        %0 of a segment was defined as the center of the histogram bin with
        %the highest absolute frequency.
        f0 = pitch(audio_data,Fs,'WindowLength',round(0.025*Fs),...
            'OverlapLength',round(0.015*Fs),'Range',[50 1000],'Method','CEP');
        min_f0=min(f0);
        max_f0=max(f0);
        mean_f0=mean(f0);
        median_f0=median(f0);
        mode_f0=mode(f0);
        var_f0=var(f0);
        %not sure how to complete the rest
        %f0 = pitch(audio_data,Fs,'WindowLength',length(audio_data),'Range',[50 1000],'Method','CEP');
        f0=-1; %note can only take window length of max 192000 which does not work for longer signals 
        
        [ZRC] = temporal_features(audio_data,Fs); 
        % ZRC: Zero Crossing Rate

        [output_spectral_features, periodogram_pks_features, ~,~,~,~,I,S] = spectral_features(audio_data,Fs);
        output_spectral_features=output_spectral_features';
        periodogram_pks_features=periodogram_pks_features';
        % spectral_features: This function computes spectral features.

        %% INPUT AND OUTPUT
        % -- Inputs --
        % x 'audio signal'
        % Fs 'sampling frequency
        % -- Outputs --
        % outputs =  spectral features including
        % meanPSD: mean frequency of power spectrum
        % stdPSD: std frequency of power spectrum
        % medPSD:  median frequency of power spectrum
        % bw: 3db bandwidth
        % p25: 1st quartile of power spectrum
        % p75: 3rd quartile of power spectrum
        % IQR: inter quartile range
        % TP: total power in 100-1000 Hz
        % p100_200: power ratio: 100-200 hz/TP
        % p200_400: power ratio: 200-400 hz/TP
        % p400_800: power ratio: 400-800 hz/TP
        % spectrum_slope2: spectrum slope
        % r_square2: R^2 statistics (linear regression for slope)
        % periodogram_pks_features: number of peaks, frequency of these peaks and 2 higher peaks frequency differences
        %I and S are intercept and slope of Slope of the regression line (SL) in db/octave

        [output_lpc,  output_lsf] = lpc_lsf_coeff(audio_data, Fs);
        % -- Inputs --
        % y: the input signal
        % -- Outputs --
        % output_mean_lpc: mean of LPC (Linear Predictive Coding) coefficients
        % output_mean_lsf: mean of LSF (Line Spectral Frequencies) coefficients
    end

    function [sdrSQI] = get_signal_to_noise_psd_score(audio_data,Fs)
        % function [sdrSQI] = get_signal_to_noise_psd_score(audio_data, Fs)
        %
        % A function to extract the sdr signal quality index from:
        % D. Springer et al., "Automated signal quality assessment of mobile
        % phone-recorded heart sound signals," JMET, In preparation, 2016
        %
        %% INPUTS:
        % audio_data: The audio data of the heart sound recording
        % Fs: the sampling frequency of the audio recording. This should be greater
        % than 2000 Hz, as this function finds the ratio of the power between
        % 20-150 Hz and the power above 600 Hz. In order to have any power above
        % 600 Hz, the sampling frequency must be greater than 1200 Hz.
        %
        %% OUTPUTS:
        % sdrSQI: Spectral distribution ratio between the expected frequencies
        % within the fundamental heart sounds and noise
        %
        %% Copyright (C) 2016  David Springer
        % dave.springer@gmail.com
        %
        % This program is free software: you can redistribute it and/or modify
        % it under the terms of the GNU General Public License as published by
        % the Free Software Foundation, either version 3 of the License, or
        % any later version.
        %
        % This program is distributed in the hope that it will be useful,
        % but WITHOUT ANY WARRANTY; without even the implied warranty of
        % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        % GNU General Public License for more details.
        %
        % You should have received a copy of the GNU General Public License
        % along with this program.  If not, see <http://www.gnu.org/licenses/>.
        if Fs < 2000
            warning('The sampling frequency in sdrSQI is low. This may result in innacutate calculation of the signal quality index');
        end
        
        
        %Get the PSD using Welch's method
        %Set window size to 2 seconds, to contain at least one heart cycle
        %Set the overlap to 50%
        %Set the number of FFT points to roughly the Fsuency, so that we get about a 1 Hz resolution (which should be acceptable, as we are looking for difference in the spectrum far larger than 1 Hz)
        
        % Find the next greatest power of 2 from the sampling frequency
        p = nextpow2(Fs*2);
        window_size = 2^p;
        
        
        [Pxx,F] = pwelch(audio_data,window_size,window_size/2,window_size,Fs);
        
        % Look for ratio of frequencies where you expect to heart heart sounds
        % (20-150 Hz) compared to those of just noise (upwards of 600 Hz).
        
        %         [~, pos1] = (min(abs(F-20)));
        %         [~, pos2] = (min(abs(F-150)));
        
        %[~, pos1] = (min(abs(F-20)));
        [~, pos2] = (min(abs(F-200)));
        
        [~, pos3] = (min(abs(F-(600))));
        
        Pxx = 10*log10(Pxx);
        
        %signal = sum(Pxx(pos1:pos2));
        signal = sum(Pxx(1:pos2));
        
        noise = sum(Pxx(pos3:end));
        
        sdrSQI = signal/noise;
        
    end

    function [bSQI_1] = get_bSQI_1(peaks_1,peaks_2,Fs)
        % [bSQI] = get_bSQI(audio_data,Fs,pi_vector, b_matrix,figures)
        %
        % A function to extract the bSQI metric, based on the agreement between two
        % peak detectors within a tolerance.
        %
        % Implemented in:
        % D. Springer et al., "Automated signal quality assessment of mobile
        % phone-recorded heart sound signals," JMET, In preparation, 2016
        %
        %% INPUTS:
        % audio_data: the PCG data
        % Fs: the sampling frequency of the audio data
        % pi_vector, b_matrix - trained parameters for the HMM-based method. Loaded
        % from the "hmm.mat" file which was trained on different, substantial,
        % database as outlined in the paper by Springer et al.
        % figures: optional boolean variable dictating the display of figures
        %
        %% OUTPUTS:
        % bSQI: the evel of agreement between two beat detectors
        %
        %% Copyright (C) 2016  David Springer
        % dave.springer@gmail.com
        %
        % This program is free software: you can redistribute it and/or modify
        % it under the terms of the GNU General Public License as published by
        % the Free Software Foundation, either version 3 of the License, or
        % any later version.
        %
        % This program is distributed in the hope that it will be useful,
        % but WITHOUT ANY WARRANTY; without even the implied warranty of
        % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        % GNU General Public License for more details.
        %
        % You should have received a copy of the GNU General Public License
        % along with this program.  If not, see <http://www.gnu.org/licenses/>.
        %% Find agreement between detectors:
        if(isempty(peaks_1) || isempty(peaks_2))
            F1_score  = 0;
        else
            [F1_score] = Bxb_compare(peaks_1, peaks_2, 0.1*Fs);
        end
        
        bSQI_1 = F1_score;
    end
    function peak_pos=heart_peaks(qt,springer_options,Fs,audio)
        homomorphic_envelope = Homomorphic_Envelope_with_Hilbert(audio, Fs);
        homomorphic_envelope = resample(homomorphic_envelope,springer_options.audio_segmentation_Fs, Fs);
        homomorphic_envelope = normalise_signal(homomorphic_envelope);
        
        %% Find mid-position of qt segments
        %Split into s1 sounds, s2 sounds and FHSounds
        %S1 sounds
        s1_segments = (qt ==1);
        s2_segments = (qt ==3);
        
        end_points_s1 = find(diff(s1_segments));
        if(mod(length(end_points_s1),2))
            end_points_s1 = [end_points_s1, length(homomorphic_envelope)];
        end
        
        
        end_points_s2 = find(diff(s2_segments));
        if(mod(length(end_points_s2),2))
            
            end_points_s2 = [end_points_s2, length(homomorphic_envelope)];
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
        peak_pos = round(((sort([mid_points_s1 mid_points_s2]))./springer_options.audio_segmentation_Fs).*Fs);
        
    end
% function SVD_SQI = get_SVD_score(truncated_autocorrelation,Fs)
%
% A function to extract the SVD signal quality index from the autocorrelation of a heart sound recording.
% This index is derived from the paper:
% D. Kumar et al., ?Noise detection during heart sound recording using
% periodicity signatures.,? Physiol. Meas., vol. 32, no. 5, pp. 599?618,
% May 2011.

% and implemented in:
% D. Springer et al., "Automated signal quality assessment of mobile
% phone-recorded heart sound signals," JMET, In preparation, 2016
%
%% INPUTS:
% truncated_autocorrelation: The truncated autocorrelation found in the
% function "get_autocorrelation.m"
% Fs: the sampling frequency of the autocorrelation
%
%% OUTPUTS:
% SVD_SQI: the ratio of second to first singular values of the SVD of
% windows of the autocorrelation
%
%% Copyright (C) 2016  David Springer
% dave.springer@gmail.com
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

    function SVD_SQI = get_SVD_score(truncated_autocorrelation,Fs,max_HR,min_HR)
        
        rho = [];
        
        
        % Set limits on the size of the window.
        % A sample of (215/500)*Fs corresponds to a heartrate of 140 BPM
        % (2)*Fs corresponds to 30 BPM:
        %min_HR=30
        %max_HR=140
        
        %max_HR=220;
        %min_HR=70;
        
        start_window_size = round((60/max_HR)*Fs);
        stop_window_size = round((60/min_HR)*Fs);
        
        % Split the autocorrelation into stacked windows
        % Find the SVD
        % Find the squared ratio of the second to first singular values:
        for T = start_window_size:5:stop_window_size
            Y = [];
            for window = 0: floor(length(truncated_autocorrelation)/T)-1
                Y = [Y;truncated_autocorrelation(window*T +1: T*(window+1))];
            end
            S_1 = svd(Y');
            if(numel(S_1) == 1)
                rho = [rho,10];
            else
                rho = [rho, ((S_1(2))/(S_1(1)))^2];
            end
            Y = [];
        end
        
        % Then, return the minimum ratio found:
        SVD_SQI = min(rho);
        
    end

% function [ccSQI] = get_ccSQI(audio_data, Fs,figures)
%
% A function to extract the cos correlation signal quality index from:
% D. Springer et al., "Automated signal quality assessment of mobile
% phone-recorded heart sound signals," JMET, In preparation, 2016
%
%% INPUTS:
% untruncated_autocorrelation: The untruncated autocorrelation of the audio
% data from the PCG recording. Derived in "get_autocorrelation.m".
% Fs: the sampling frequency of the audio recording
% figures: (optional) boolean variable for displaying figures
%
%% OUTPUTS:
% ccSQI: Correlation between the autocorrelation signal and a fitted sinusoid
%
%% Copyright (C) 2016  David Springer
% dave.springer@gmail.com
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


    function [ccSQI] = get_ccSQI(untruncated_autocorrelation, Fs, max_HR,min_HR,figures)
        
        if(nargin < 5)
            figures = false;
        end
        
        if(length(untruncated_autocorrelation) < 5*Fs)
            ccSQI = 0;
            warning('PCG signal too short - cannot comput ccSQI');
        else
            
            t = 0:1:length(untruncated_autocorrelation)-1; % Time Samples
            sums = [];
            
            % potential range of heart rates (in cycles per second rather than
            % BPM):
            %min_HR=0.6*60; %36
            %max_HR=2.33*60; %139.8
            %min_HR=70;
            %max_HR=220;
            frequency_range = (min_HR/60):0.01:(max_HR/60);
            
            for i = 1:length(frequency_range)
                
                f = frequency_range(i);% Heart beat frequency
                
                % Find the samples each side of the cosine peaks. These peaks are
                % at multiples of the "f" frequency of the cosine:
                cos_peaks1 = (1:5)'.*Fs/f - round(Fs*0.12);
                cos_peaks2 = (1:5)'.*Fs/f + round(Fs*0.12);
                cos_spread = [];
                for j = 1:5
                    % Find the "spread" of the samples each side of each cosine
                    % peak - cos_peaks1 is 0.12*Fs samples to the left of the peak
                    % cos_peaks2 is 0.12*Fs to the right of each peak, and the
                    % cos_spread finds the samples in between these two points, for
                    % the five cosine peaks. This is unless the spread passes the
                    % length of the signal. In that case, the spread is truncated
                    % to the length.
                    cos_spread = [cos_spread, round(min([cos_peaks1(j) length(untruncated_autocorrelation)]):min([cos_peaks2(j) length(untruncated_autocorrelation)]))];
                end
                % Find the squared sum of the samples of the autocorrelation
                % function within the "cos_spread" each side of the cosine peaks
                sums = [sums sum((untruncated_autocorrelation(cos_spread)).^2)];
            end
            
            % Find the heart frequency which resulted in the maximum sum of the
            % autocorrelation samples within the "cos_spread" of the fitted cosine:
            [~,b] = max(sums);
            f = frequency_range(b);
            oscil = (cos(2*pi*f/Fs*t));
            oscil = ((oscil)>0).*(oscil);
            
            % Limit the correlation to 5 seconds:
            max_cos_peaks = round(5*Fs/f) + round(Fs*0.120);
            correlation_coefficient = corrcoef(untruncated_autocorrelation(1:max_cos_peaks), oscil(1:max_cos_peaks));
            
            ccSQI =  correlation_coefficient(2,1);
            
            if(figures)
                figure('Name', 'ccSQI Plot');
                t = (1:length(untruncated_autocorrelation(1:max_cos_peaks)))./Fs;
                plot(t,untruncated_autocorrelation(1:max_cos_peaks));
                hold on;
                plot(t,oscil(1:max_cos_peaks),'r');
            end
        end
        
        
        
    end

% function [mSQI] = get_peak_in_autocorrelation(untruncated_autocorrelation, Fs, figures)
%
% A function to extract amplitude of the maximum peak in the autocorrelation
% signal between lag values of 0.5 and 2.33 seconds.
% Developed as part of:
% D. Springer et al., "Automated signal quality assessment of mobile
% phone-recorded heart sound signals," JMET, In preparation, 2016
%
% This function requires the Matlab implementation of the "sampen.m" function from physionet:
% https://physionet.org/physiotools/sampen/
%
%% INPUTS:
% untruncated_autocorrelation: The untruncated autocorrelation of the audio
% data from the PCG recording. Derived in "get_autocorrelation.m".
% Fs: the sampling frequency of the autocorrelation
% figures: optional boolean variable to display figures
%
%% OUTPUTS:
% sample_entropy: Sample entropy of the autocorrelation signal
%
%% Copyright (C) 2016  David Springer
% dave.springer@gmail.com
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

    function [mSQI, pos,corr_prominance] = get_max_peak_in_autocorrelation(untruncated_autocorrelation,Fs,max_HR,min_HR, figures)
        if nargin < 5
            figures = 0;
        end
        
        % Set the heart rate limits
        % 0.43*Fs corresponds to 140 BPM
        % 2*Fs corresponds to 30 BPM
        %min_HR=30;
        %max_HR=140;
        %min_HR=70;
        %max_HR=220;
        
        start_point = round(Fs*60/max_HR);
        end_point = round(Fs*60/min_HR);
        
        % Find the amplitude and position of the maximum peak between these two limits:
        [mSQI, pos] = max(untruncated_autocorrelation(start_point:end_point));
        pos1=pos;
        pos = pos + start_point + 1;
        
        %correlation prominance
        % prominance value of the max peak found in mSQI
        [~,locs,~,p] = findpeaks(untruncated_autocorrelation(start_point:end_point));
        corr_prominance=p(locs==pos1);
        if isempty(corr_prominance)
            corr_prominance=1;
        end
        
        
        if(figures)
            figure('Name','mSQI plot');
            plot(untruncated_autocorrelation);
            hold on;
            plot(pos,mSQI,'ko');
            legend('Autocorrelation', 'Position of maximum peak between heart rate limits');
        end
        
        
    end
end