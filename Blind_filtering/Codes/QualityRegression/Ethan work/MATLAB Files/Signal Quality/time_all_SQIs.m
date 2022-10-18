function  t= time_all_SQIs(audio_data, Fs)
SQI=table;
max_HR=220;
min_HR=70;
max_RR=80;
min_RR=15;

%dealing with padded zeros and instances of zeros
audio_data(audio_data==0)=min(abs(audio_data(audio_data~=0))); 
audio_data=resample(audio_data,4000,Fs);
f = @() resample(audio_data,4000,Fs);
t.pre_4000=timeit(f);
Fs=4000;

Fs_2000=2000; 
audio_downsampled_2000 = resample(audio_data,Fs_2000,Fs);

Fs_1000=1000;
audio_downsampled_1000 = resample(audio_data,Fs_1000,Fs);

%% Audio derived features %%
% Sample Entropy-----------------------------------------------------------
% Improving the quality of point of care diagnostics with real-time machine learning in low literacy LMIC setting
% https://dl.acm.org/doi/pdf/10.1145/3209811.3209815
%[SQI.audio_se] = get_entropy(audio_data, 2, 0.1);
%f = @() get_entropy(audio_data, 2, 0.1);
%t.audio_se=timeit(f);
audio_downsampled_30=resample(audio_data,30,Fs);
%SQI.audio_se = get_entropy(audio_downsampled_30, 2, 0.1);

f = @() resample(audio_data,30,Fs);
t.audio_se1=timeit(f);
f = @() get_entropy(audio_downsampled_30, 2, 0.1);
t.audio_se2=timeit(f);

% Clipping-----------------------------------------------------------------
% Computerized lung sound screening for pediatric auscultation in noisy field environments
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7953509
%SQI.audio_clip=detect_clipping(audio_data);
f = @() detect_clipping(audio_data);
t.audio_clip=timeit(f);

% Rate average energy E[R(f)]----------------------------------------------
% An objective measure of signal quality for pediatric lung auscultations
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9176539
%SQI.audio_rate_average_energy_2000=get_rate_average_energy(audio_downsampled_2000, Fs_2000); 
%SQI.audio_rate_average_energy_1000=get_rate_average_energy(audio_downsampled_1000, Fs_1000); 
f = @() get_rate_average_energy(audio_downsampled_2000, Fs_2000); 
t.audio_rate_average_energy_2000=timeit(f);
f = @() get_rate_average_energy(audio_downsampled_1000, Fs_1000); 
t.audio_rate_average_energy_1000=timeit(f);

% Heart contamination------------------------------------------------------
% Computerized lung sound screening for pediatric auscultation in noisy field environments
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7953509
%[SQI.audio_heart_con1,SQI.audio_heart_con2]=get_heart_contamination(audio_data,Fs); 
f = @() get_heart_contamination(audio_data,Fs); 
t.audio_heart_con=timeit(f);

% Skewness-----------------------------------------------------------------
% Automatic signal quality index determination of radar-recorded heart sound signals using ensemble classification
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8731709
%[SQI.audio_skew] = get_skewness_score(audio_data);
f = @() get_skewness_score(audio_data);
t.audio_skew=timeit(f);


% Kurtosis-----------------------------------------------------------------
% Automated signal quality assessment of mobile phone-recorded heart sound signals
% https://www.tandfonline.com/doi/pdf/10.1080/03091902.2016.1213902?needAccess=true;
%[SQI.audio_kurtosis] = get_kurtosis_score(audio_data);
f = @() get_kurtosis_score(audio_data);
t.audio_kurtosis=timeit(f);

% Variance-----------------------------------------------------------------
% Improving the quality of point of care diagnostics with real-time machine learning in low literacy LMIC setting
% https://dl.acm.org/doi/pdf/10.1145/3209811.3209815
%SQI.audio_var = get_variance_score(audio_data);
f = @() get_variance_score(audio_data);
t.audio_var=timeit(f);

% High frequency variance--------------------------------------------------
% Detection of pathological heart sounds
% https://iopscience.iop.org/article/10.1088/1361-6579/aa7840/pdf
% Computerized quality assessment of phonocardiogram signal measurement-acquisition parameters
% https://www.researchgate.net/publication/225085127_Computerized_quality_assessment_of_phonocardiogram_signal_measurement-acquisition_parameters
audio_high =butterworth_high_pass_filter(audio_data,2,700,Fs);
%SQI.audio_high_var = get_variance_score(audio_high);
f = @() butterworth_high_pass_filter(audio_data,2,700,Fs);
t.audio_high=timeit(f);
f = @() get_variance_score(audio_high);
t.audio_high_var=timeit(f);

% Zero crossing rate-------------------------------------------------------
% Preterm Neonates Breathing Sounds Analysis
% https://github.com/JulieKy/Internship-2A
%SQI.audio_zcr = get_zrc(audio_data);
f = @() get_zrc(audio_data);
t.audio_zcr=timeit(f);


% Linear predictive coefficient--------------------------------------------
% Heart sound anomaly and quality detection using ensemble of neural networks without segmentation
% https://ieeexplore.ieee.org/document/7868817
% Linear Predictive Coefficient (LPC): 10th order linear predictor with all coefficients are used as features
%[SQI.audio_lpc10,~] = lpc(double(audio_data),10);
f = @() lpc(double(audio_data),10);
t.audio_lpc10=timeit(f);

% Preterm Neonates Breathing Sounds Analysis
% https://github.com/JulieKy/Internship-2A
%[SQI.audio_lpc6,  SQI.audio_lsf6] = lpc_lsf_coeff_modified(audio_data);
f = @() lpc_lsf_coeff_modified(audio_data);
t.audio_lpc_lsf_10=timeit(f);

% Entropy based features---------------------------------------------------
% Heart sound anomaly and quality detection using ensemble of neural networks without segmentation
% https://ieeexplore.ieee.org/document/7868817
%[SQI.audio_shannon,SQI.audio_tsallis,SQI.audio_renyi]=get_prob_entropy(audio_data);
f = @() get_prob_entropy(audio_data);
t.audio_prob_entropy=timeit(f);

% Degree of periodicity----------------------------------------------------
% Automated Signal Quality Assessment for Heart Sound Signal by Novel Features and Evaluation in Open Public Datasets
% https://downloads.hindawi.com/journals/bmri/2021/7565398.pdf
% Best subsequence selection of heart sound recording based on degree of sound periodicity
% https://www.researchgate.net/publication/224249101_Best_subsequence_selection_of_heart_sound_recording_based_on_degree_of_sound_periodicity
% cycle frequency considered
%Min_cf=0.3;
%Max_cf=2.5;
Min_cf = (60/max_HR);
Max_cf = (60/min_HR);
%[SQI.audio_periodicity_heart,~,~]=getDegree_cycle_modified(audio_data,Min_cf,Max_cf,Fs);
f = @() getDegree_cycle_modified(audio_data,Min_cf,Max_cf,Fs);
t.audio_periodicity_heart=timeit(f);

Min_cf = (60/max_RR);
Max_cf = (60/min_RR);
%[SQI.audio_periodicity_lung,~,~]=getDegree_cycle_modified(audio_data,Min_cf,Max_cf,Fs);
f = @() getDegree_cycle_modified(audio_data,Min_cf,Max_cf,Fs);
t.audio_periodicity_lung=timeit(f);

% Kurtosis-----------------------------------------------------------------
% Improving the quality of point of care diagnostics with real-time machine learning in low literacy LMIC setting
% https://dl.acm.org/doi/pdf/10.1145/3209811.3209815
% Find the autocorrelation:
acf = autocorr(audio_data,length(audio_data)-1, [] , 2);
%[SQI.audio_autocorr_kurtosis] = get_kurtosis_score(acf);
f = @()  autocorr(audio_data,length(audio_data)-1, [] , 2);
t.audio_autocorr=timeit(f);
f = @() get_kurtosis_score(acf);
t.audio_autocorr_kurtosis=timeit(f);

% Wavelet features---------------------------------------------------------
% Improving the quality of point of care diagnostics with real-time machine learning in low literacy LMIC setting
% https://dl.acm.org/doi/pdf/10.1145/3209811.3209815
%[SQI.wav_energy] = get_wavelet_features(audio_data);
f = @() get_wavelet_features(audio_data);
t.wav_energy=timeit(f);

% Heart sound anomaly and quality detection using ensemble of neural networks without segmentation
% https://ieeexplore.ieee.org/document/7868817
%[SQI.wav_log, SQI.wav_shannon_a5,SQI.wav_tsallis_a5,SQI.wav_renyi_a5,...
%    SQI.wav_shannon_d5,SQI.wav_tsallis_d5, SQI.wav_renyi_d5,SQI.wav_shannon_d4,...
%    SQI.wav_tsallis_d4,SQI.wav_renyi_d4]=get_wavelet_features_zabihi(audio_data);
f = @() get_wavelet_features_zabihi(audio_data);
t.wav_zabihi=timeit(f);

% Analysis of PCG signals using quality assessment and homomorphic filters for localization and classification of heart sounds
% https://reader.elsevier.com/reader/sd/pii/S0169260718303596?token=BCE2CE37A7B30CA802618D72F26C8110294D3BE20C6C531A50359E2E61DCE98550072BF9C21048BBCF075141E50D1C09&originRegion=us-east-1&originCreation=20210526053544
%[SQI.wav_RMSSD,  SQI.wav_RZC]=get_wavelet_mubarak(audio_data);
f = @() get_wavelet_mubarak(audio_data);
t.wav_mubarak=timeit(f);

% Fundamental frequency----------------------------------------------------
% Improving the quality of point of care diagnostics with real-time machine learning in low literacy LMIC setting
% https://dl.acm.org/doi/pdf/10.1145/3209811.3209815
%[SQI.f0_overall,SQI.f0_min,SQI.f0_max,SQI.f0_mean,SQI.f0_median,SQI.f0_mode,SQI.f0_var,SQI.cry_f0]...
%    = get_fundemental_frequency(audio_data,Fs);
f = @() get_fundemental_frequency(audio_data,Fs);
t.f0=timeit(f);

%% Get Autocorrelation functions
% Remove noise spikes in the signal using the method developed by Schmidt et al:
% Segmentation of heart sound recordings by a duration-dependent hidden Markov model
% https://iopscience.iop.org/article/10.1088/0967-3334/31/4/004/pdf?casa_token=77qaFGzyglgAAAAA:-Iz5cUPtusnih53uEpJOySYpDfpE5WGbDzUrAJFApWe7MuTI9LEfNUHA1gMacLQSkF3k2ayY2HM
audio_downsampled_2000 = schmidt_spike_removal(audio_downsampled_2000,Fs_2000);
[truncated_autocorrelation, untruncated_autocorrelation, semi_truncated_autocorrelation] = get_autocorrelation_modified(audio_downsampled_2000, Fs_2000);

f = @() schmidt_spike_removal(audio_downsampled_2000,Fs_2000);
t.audio_downsampled=timeit(f);
f = @() get_autocorrelation_modified(audio_downsampled_2000, Fs_2000);
t.autocorr=timeit(f);
%% Autocorrelation derived features
% Sample Entropy-----------------------------------------------------------
% Automated signal quality assessment of mobile phone-recorded heart sound signals
% https://www.tandfonline.com/doi/pdf/10.1080/03091902.2016.1213902?needAccess=true;
%SQI.trunc_autocorr_se1 = get_entropy(truncated_autocorrelation, 2, 0.0008);
%SQI.strunc_autocorr_se1 = get_entropy(truncated_autocorrelation, 2, 0.0008);
f = @() get_entropy(truncated_autocorrelation, 2, 0.0008);
t.trunc_autocorr_se1=timeit(f);
f = @() get_entropy(truncated_autocorrelation, 2, 0.0008);
t.strunc_autocorr_se1=timeit(f);

% Signal quality classification of mobile phone-recorded phonocardiogram signals
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6853814
%M=1,r=0.01
%SQI.trunc_autocorr_se2 = get_entropy(truncated_autocorrelation, 1, 0.01);
%SQI.strunc_autocorr_se2 = get_entropy(semi_truncated_autocorrelation, 1, 0.01);
f = @() get_entropy(truncated_autocorrelation, 1, 0.01);
t.trunc_autocorr_se2=timeit(f);
f = @() get_entropy(truncated_autocorrelation, 1, 0.01);
t.strunc_autocorr_se2=timeit(f);
%M=2, r=0.01
%SQI.trunc_autocorr_se3 = get_entropy(truncated_autocorrelation, 2, 0.01);
%SQI.strunc_autocorr_se3 = get_entropy(semi_truncated_autocorrelation, 2, 0.01);
f = @() get_entropy(truncated_autocorrelation, 2, 0.01);
t.trunc_autocorr_se3=timeit(f);
f = @() get_entropy(truncated_autocorrelation, 2, 0.01);
t.strunc_autocorr_se3=timeit(f);
%M=1, r=0.001
%SQI.trunc_autocorr_se4 = get_entropy(truncated_autocorrelation, 1, 0.001);
%SQI.strunc_autocorr_se4 = get_entropy(semi_truncated_autocorrelation, 1, 0.001);
f = @() get_entropy(truncated_autocorrelation, 1, 0.001);
t.trunc_autocorr_se4=timeit(f);
f = @() get_entropy(truncated_autocorrelation, 1, 0.001);
t.strunc_autocorr_se4=timeit(f);
%M=2, r=0.001
%SQI.trunc_autocorr_se5 = get_entropy(truncated_autocorrelation, 2, 0.001);
%SQI.strunc_autocorr_se5 = get_entropy(semi_truncated_autocorrelation, 2, 0.001);
f = @() get_entropy(truncated_autocorrelation, 2, 0.001);
t.trunc_autocorr_se5=timeit(f);
f = @() get_entropy(truncated_autocorrelation, 2, 0.001);
t.strunc_autocorr_se5=timeit(f);
% Automated Signal Quality Assessment for Heart Sound Signal by Novel Features and Evaluation in Open Public Datasets
% https://downloads.hindawi.com/journals/bmri/2021/7565398.pdf
% M=2, r=0.2
truncated_autocorrelation_down=resample(truncated_autocorrelation,30,Fs_2000);
semi_truncated_autocorrelation_down=resample(semi_truncated_autocorrelation,30,Fs_2000);
%SQI.trunc_autocorr_se6 = get_entropy(truncated_autocorrelation_down, 2, 0.2);
%SQI.strunc_autocorr_se6 = get_entropy(semi_truncated_autocorrelation_down, 2, 0.2);
f = @() resample(truncated_autocorrelation,30,Fs_2000);
t.trunc_autocorr_30=timeit(f);
f = @() resample(semi_truncated_autocorrelation,30,Fs_2000);
t.strunc_autocorr_30=timeit(f);
f = @() get_entropy(truncated_autocorrelation_down, 2, 0.2);
t.trunc_autocorr_se6=timeit(f);
f = @() get_entropy(semi_truncated_autocorrelation_down, 2, 0.2);
t.strunc_autocorr_se6=timeit(f);

% Hjorth activity score----------------------------------------------------
% Automated signal quality assessment of mobile phone-recorded heart sound signals
% https://www.tandfonline.com/doi/pdf/10.1080/03091902.2016.1213902?needAccess=true;
%SQI.trunc_autocorr_hjorth = get_hjorth_activity_score(truncated_autocorrelation);
%SQI.strunc_autocorr_hjorth = get_hjorth_activity_score(semi_truncated_autocorrelation);
f = @() get_hjorth_activity_score(truncated_autocorrelation);
t.trunc_autocorr_hjorth=timeit(f);
f = @() get_hjorth_activity_score(semi_truncated_autocorrelation);
t.strunc_autocorr_hjorth=timeit(f);

% Variance-----------------------------------------------------------------
% Automated signal quality assessment of mobile phone-recorded heart sound signals
% https://www.tandfonline.com/doi/pdf/10.1080/03091902.2016.1213902?needAccess=true;
%SQI.trunc_autocorr_var = get_variance_score(truncated_autocorrelation);
%SQI.strunc_autocorr_var = get_variance_score(semi_truncated_autocorrelation);
f = @() get_variance_score(truncated_autocorrelation);
t.trunc_autocorr_var=timeit(f);
f = @() get_variance_score(semi_truncated_autocorrelation);
t.strunc_autocorr_var=timeit(f);

% Singular value decomposition---------------------------------------------
% Automated signal quality assessment of mobile phone-recorded heart sound signals
% https://www.tandfonline.com/doi/pdf/10.1080/03091902.2016.1213902?needAccess=true;
% Noise detection during heart sound recording using periodicity signatures
% https://www.researchgate.net/publication/51037539_Noise_detection_during_heart_sound_recording_using_periodicity_signatures
%SQI.trunc_autocorr_svd_heart = get_SVD_score_modified(truncated_autocorrelation,max_HR,min_HR,Fs_2000);
%SQI.strunc_autocorr_svd_heart = get_SVD_score_modified(semi_truncated_autocorrelation,max_HR,min_HR,Fs_2000);
%SQI.trunc_autocorr_svd_lung = get_SVD_score_modified(truncated_autocorrelation,max_RR,min_RR,Fs_2000);
%SQI.strunc_autocorr_svd_lung = get_SVD_score_modified(semi_truncated_autocorrelation,max_RR,min_RR,Fs_2000);
f = @() get_SVD_score_modified(truncated_autocorrelation,max_HR,min_HR,Fs_2000);
t.trunc_autocorr_svd_heart =timeit(f);
f = @() get_SVD_score_modified(semi_truncated_autocorrelation,max_HR,min_HR,Fs_2000);
t.strunc_autocorr_svd_heart=timeit(f);
f = @() get_SVD_score_modified(truncated_autocorrelation,max_RR,min_RR,Fs_2000);
t.trunc_autocorr_svd_lung=timeit(f);
f = @() get_SVD_score_modified(semi_truncated_autocorrelation,max_RR,min_RR,Fs_2000);
t.strunc_autocorr_svd_lung=timeit(f);

% Correlation with cosine--------------------------------------------------
% Automated signal quality assessment of mobile phone-recorded heart sound signals
% https://www.tandfonline.com/doi/pdf/10.1080/03091902.2016.1213902?needAccess=true;
%SQI.autocorr_cc_heart = get_ccSQI_modified(untruncated_autocorrelation,max_HR,min_HR,Fs_2000);
%SQI.autocorr_cc_lung = get_ccSQI_modified(untruncated_autocorrelation,max_RR,min_RR,Fs_2000);
f = @() get_ccSQI_modified(untruncated_autocorrelation,max_HR,min_HR,Fs_2000);
t.autocorr_cc_heart=timeit(f);
f = @() get_ccSQI_modified(untruncated_autocorrelation,max_RR,min_RR,Fs_2000);
t.autocorr_cc_lung=timeit(f);

% Max autocorrelation------------------------------------------------------
% Automated signal quality assessment of mobile phone-recorded heart sound signals
% https://www.tandfonline.com/doi/pdf/10.1080/03091902.2016.1213902?needAccess=true;
% Automatic signal quality index determination of radar-recorded heart sound signals using ensemble classification
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8731709
% Automated Signal Quality Assessment for Heart Sound Signal by Novel Features and Evaluation in Open Public Datasets
% https://downloads.hindawi.com/journals/bmri/2021/7565398.pdf
%[SQI.autocorr_max_val_heart, SQI.autocorr_cycle_dur_heart, SQI.autocorr_max_prom_heart] =...
%    get_max_peak_in_autocorrelation_modified(untruncated_autocorrelation,max_HR,min_HR,Fs_2000);
%[SQI.autocorr_max_val_lung, SQI.autocorr_cycle_dur_lung, SQI.autocorr_max_prom_lung] =...
%    get_max_peak_in_autocorrelation_modified(untruncated_autocorrelation,max_RR,min_RR,Fs_2000);
f = @() get_max_peak_in_autocorrelation_modified(untruncated_autocorrelation,max_HR,min_HR,Fs_2000);
t.autocorr_max_heart=timeit(f);
f = @() get_max_peak_in_autocorrelation_modified(untruncated_autocorrelation,max_RR,min_RR,Fs_2000);
t.autocorr_max_lung=timeit(f);

% Kurtosis-----------------------------------------------------------------
% Improving the quality of point of care diagnostics with real-time machine learning in low literacy LMIC setting
% https://dl.acm.org/doi/pdf/10.1145/3209811.3209815
%[SQI.autocorr_kurtosis] = get_kurtosis_score(untruncated_autocorrelation);
f = @() get_kurtosis_score(untruncated_autocorrelation);
t.autocorr_kurtosis=timeit(f);

%% Get Power Spectral Density
% Automated signal quality assessment of mobile phone-recorded heart sound signals
% https://www.tandfonline.com/doi/pdf/10.1080/03091902.2016.1213902?needAccess=true;
[Pxx,F] = get_psd(audio_downsampled_2000,Fs_2000);
f = @() get_psd(audio_downsampled_2000,Fs_2000);
t.power=timeit(f);
%% Power spectral density derived features
% Peaks features-----------------------------------------------------------
% Preterm Neonates Breathing Sounds Analysis
% https://github.com/JulieKy/Internship-2A
pxx_smooth=get_maf(Pxx);
nb_higherPks_MAF=2;
%[SQI.power_maf_num_peaks,  SQI.power_maf_f_peaks, SQI.power_maf_diff_peaks]...
%    = peaks_features_modified(pxx_smooth,F, nb_higherPks_MAF);
f = @() get_maf(Pxx);
t.power_maf1=timeit(f);
f = @() peaks_features_modified(pxx_smooth,F, nb_higherPks_MAF);
t.power_maf2=timeit(f);

[fi_tot,SQI.power_gmm_parameters]=get_gmm(Pxx,F);
nb_higherPks_GMM=2;
%[SQI.power_gmm_num_peaks,  SQI.power_gmm_f_peaks, SQI.power_gmm_diff_peaks]...
%    = peaks_features_modified(fi_tot,F, nb_higherPks_GMM);
f = @() get_gmm(Pxx,F);
t.power_gmm1=timeit(f);
f = @() peaks_features_modified(fi_tot,F, nb_higherPks_GMM);
t.power_gmm2=timeit(f);

% Linear regression line--------------------------------------------------- 
% Preterm Neonates Breathing Sounds Analysis
% https://github.com/JulieKy/Internship-2A
%[SQI.power_regression_slope,SQI.power_regression_intercept,SQI.power_regression_r2]=get_regression_line(Pxx,F);
f = @() get_regression_line(Pxx,F);
t.power_regression=timeit(f);

% Signal to noise ratio----------------------------------------------------
% Automated signal quality assessment of mobile phone-recorded heart sound signals
% https://www.tandfonline.com/doi/pdf/10.1080/03091902.2016.1213902?needAccess=true;
%SQI.power_sdr1_heart = get_signal_to_noise_psd_score_modified(Pxx,F,20,150,600);
%SQI.power_sdr2_heart = get_signal_to_noise_psd_score_modified(Pxx,F,20,200,600);
%SQI.power_sdr1_lung = get_signal_to_noise_psd_score_modified(Pxx,F,100,700,1000);
%SQI.power_sdr2_lung = get_signal_to_noise_psd_score_modified(Pxx,F,100,700,0);
f = @() get_signal_to_noise_psd_score_modified(Pxx,F,20,150,600);
t.power_sdr1_heart=timeit(f);
f = @() get_signal_to_noise_psd_score_modified(Pxx,F,20,200,600);
t.power_sdr2_heart=timeit(f);
f = @() get_signal_to_noise_psd_score_modified(Pxx,F,100,700,1000);
t.power_sdr1_lung=timeit(f);
f = @() get_signal_to_noise_psd_score_modified(Pxx,F,100,700,0);
t.power_sdr2_lung=timeit(f);

% Computerized quality assessment of phonocardiogram signal measurement-acquisition parameters
% https://www.researchgate.net/publication/225085127_Computerized_quality_assessment_of_phonocardiogram_signal_measurement-acquisition_parameters
%SQI.power_sdr3_heart = get_signal_to_noise_psd_score_modified(Pxx,F,20,200,0);
f = @() get_signal_to_noise_psd_score_modified(Pxx,F,20,200,0);
t.power_sdr3_heart=timeit(f);

% Cry detection------------------------------------------------------------
% Preterm Neonates Breathing Sounds Analysis
% https://github.com/JulieKy/Internship-2A
%SQI.power_sdr_cry = get_signal_to_noise_psd_score_modified(Pxx,F,295,406,0);
f = @() get_signal_to_noise_psd_score_modified(Pxx,F,295,406,0);
t.power_sdr_cry1=timeit(f);
% Computerized lung sound screening for pediatric auscultation in noisy field environments
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7953509
f = @() get_signal_to_noise_psd_score_modified(Pxx,F,1000,F(end),0);
t.power_sdr_cry2=timeit(f);

% Spectral composition-----------------------------------------------------
% Heart sound anomaly and quality detection using ensemble of neural networks without segmentation
% https://ieeexplore.ieee.org/document/7868817
%[SQI.power_b1, SQI.power_b2, SQI.power_b3, SQI.power_b4, SQI.power_b5, ...
%    SQI.power_b6, SQI.power_b7, SQI.power_b8, SQI.power_b9, SQI.power_b10] = get_spectral_composition(Pxx,F);
f = @() get_spectral_composition(Pxx,F);
t.power_b=timeit(f);

% Preterm Neonates Breathing Sounds Analysis
% https://github.com/JulieKy/Internship-2A
%[SQI.power_f_mean,SQI.power_std,SQI.power_f_med,SQI.power_bw,SQI.power_f_p25,...
%    SQI.power_f_p75,SQI.power_f_IQR,SQI.power_000_100,SQI.power_100_200,...
%    SQI.power_200_400,SQI.power_400_600,SQI.power_600_800,SQI.power_800_1000,...
%    SQI.power_1000_1200,SQI.power_1200_1400,SQI.power_1400_1600,SQI.power_1600_1800,...
%    SQI.power_1800_2000,SQI.power_tp]=get_spectrum_parameters(Pxx,F);   

f = @() get_spectrum_parameters(Pxx,F);  
t.power_spectrum=timeit(f);
        
% Power ratios-------------------------------------------------------------
% Automated Signal Quality Assessment for Heart Sound Signal by Novel Features and Evaluation in Open Public Datasets
% https://downloads.hindawi.com/journals/bmri/2021/7565398.pdf
%SQI.power_ratio_low = get_signal_to_noise_psd_score_modified(Pxx,F,24,144,0);
%SQI.power_ratio_high = get_signal_to_noise_psd_score_modified(Pxx,F,200,1000,0);
%SQI.power_ratio_higher = get_signal_to_noise_psd_score_modified(Pxx,F,1000,Fs_2000/2,0);
%SQI.power_ratio_mid = get_signal_to_noise_psd_score_modified(Pxx,F,144,200,0);

f = @() get_signal_to_noise_psd_score_modified(Pxx,F,24,144,0); 
t.power_ratio_low=timeit(f);
f = @() get_signal_to_noise_psd_score_modified(Pxx,F,200,1000,0); 
t.power_ratio_high=timeit(f);
f = @() get_signal_to_noise_psd_score_modified(Pxx,F,1000,Fs_2000/2,0);  
t.power_ratio_higher=timeit(f);
f = @() get_signal_to_noise_psd_score_modified(Pxx,F,144,200,0);
t.power_ratio_mid=timeit(f);


% Modified power spectral density centroid.--------------------------------
% Heart sound anomaly and quality detection using ensemble of neural networks without segmentation
% https://ieeexplore.ieee.org/document/7868817
%SQI.power_freq_centroid = (sum(F'.*(Pxx.^2)))/(sum(Pxx.^2));

f = @() (sum(F'.*(Pxx.^2)))/(sum(Pxx.^2));
t.power_freq_centroid=timeit(f);

% Spectral entropy---------------------------------------------------------
% Automatic signal quality index determination of radar-recorded heart sound signals using ensemble classification
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8731709
%SQI.power_hn = get_spectral_entropy(Pxx);

f = @() get_spectral_entropy(Pxx);
t.power_hn=timeit(f);

% Dominant frequency freatures---------------------------------------------
% Automatic signal quality index determination of radar-recorded heart sound signals using ensemble classification
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8731709
%[SQI.power_max_pow, SQI.power_max_freq,SQI.power_ratio_max] = get_dominant_frequency_features(Pxx,F);

f = @() get_dominant_frequency_features(Pxx,F);
t.power_max=timeit(f);

% Variance-----------------------------------------------------------------
% Improving the quality of point of care diagnostics with real-time machine learning in low literacy LMIC setting
% https://dl.acm.org/doi/pdf/10.1145/3209811.3209815
%SQI.power_var = get_variance_score(Pxx);

f = @() get_variance_score(Pxx);
t.power_var=timeit(f);

%% Get Segmentations
% Heart sound segmentation algorithm based on heart sound envelogram
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=647841
[peak_pos_liang,~,~] = get_Liang_peaks_modified(audio_downsampled_1000,Fs_1000);

f = @() get_Liang_peaks_modified(audio_downsampled_1000,Fs);
t.seg_liang=timeit(f);

% Segmentation of heart sound recordings by a duration-dependent hidden Markov model
% https://iopscience.iop.org/article/10.1088/0967-3334/31/4/004/pdf?casa_token=77qaFGzyglgAAAAA:-Iz5cUPtusnih53uEpJOySYpDfpE5WGbDzUrAJFApWe7MuTI9LEfNUHA1gMacLQSkF3k2ayY2HM
[heartRate, systolicTimeInterval] = getHeartRateSchmidt_modified(audio_downsampled_1000,max_HR,min_HR, Fs_1000);
[peak_pos_schmidt,~] = get_schmidt_hmm_peakpos_modified(audio_downsampled_1000,heartRate, systolicTimeInterval,Fs_1000);

f = @()  getHeartRateSchmidt_modified(audio_downsampled_1000,max_HR,min_HR, Fs_1000);
t.seg_heartRate=timeit(f);
f = @() get_schmidt_hmm_peakpos_modified(audio_downsampled_1000,heartRate, systolicTimeInterval,Fs_1000);
t.seg_schmidt=timeit(f);

% Logistic regression-HSMM-based heart sound segmentation
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7234876
[peak_pos_springer,assigned_states] = get_springer_hmm_peakpos_modified(audio_downsampled_1000, Fs_1000,heartRate, systolicTimeInterval);

f = @() get_springer_hmm_peakpos_modified(audio_downsampled_1000, Fs_1000,heartRate, systolicTimeInterval);
t.seg_springer=timeit(f);

% 25-400Hz 4th order Butterworth band pass
% Logistic regression-HSMM-based heart sound segmentation
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7234876
denoised_signal = butterworth_low_pass_filter(audio_downsampled_1000,2,400,Fs_1000, false);
denoised_signal = butterworth_high_pass_filter(denoised_signal,2,25,Fs,false);

f = @() butterworth_low_pass_filter(audio_downsampled_1000,2,400,Fs_1000, false);
t.seg_denoised1=timeit(f);
f = @() butterworth_high_pass_filter(denoised_signal,2,25,Fs,false);
t.seg_denoised2=timeit(f);

% Spike removal from the original paper:
% Segmentation of heart sound recordings by a duration-dependent hidden Markov model
% https://iopscience.iop.org/article/10.1088/0967-3334/31/4/004/pdf?casa_token=77qaFGzyglgAAAAA:-Iz5cUPtusnih53uEpJOySYpDfpE5WGbDzUrAJFApWe7MuTI9LEfNUHA1gMacLQSkF3k2ayY2HM
denoised_signal = schmidt_spike_removal(denoised_signal,Fs_1000);

f = @() schmidt_spike_removal(denoised_signal,Fs_1000);
t.seg_denoised3=timeit(f);

[signal_s1,signal_s2, signal_systole,signal_diastole]=seperate_states_2(assigned_states,denoised_signal,false);
%[SQI.percentage_bad_cycles_sys,SQI.percentage_bad_cycles_dia,...
%    SQI.percentage_bad_cycles_s1,SQI.percentage_bad_cycles_s2,SQI.percentage_bad_cycles_overall]...
%    =get_bad_cycles(assigned_states,signal_systole,signal_diastole,signal_s1,signal_s2);
f = @() seperate_states_2(assigned_states,denoised_signal,false);
t.seg_separate_states=timeit(f);
f = @() get_bad_cycles(assigned_states,signal_systole,signal_diastole,signal_s1,signal_s2);
t.seg_bad_cycles=timeit(f);

% Heart sounds quality analysis for automatic cardiac biometry applications
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=5386481
%[SQI.seg_qual,SQI.seg_qual_avg,SQI.seg_qual_inv,SQI.seg_qual_inv_avg] = get_segmentation_quality(assigned_states,signal_s1,signal_s2,Fs_1000);

f = @() get_segmentation_quality(assigned_states,signal_s1,signal_s2,Fs_1000);
t.seg_qual=timeit(f);
%% Get envelopes
% Logistic regression-HSMM-based heart sound segmentation
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7234876
[psd_envelope,psd_original_envelope,psd_full,psd_full_freq] = get_psd_envelope(denoised_signal, Fs_1000);
wavelet_envelope = get_wavelet_envelope(denoised_signal, Fs_1000);
homomorphic_envelope = get_homomorphic_envelope(denoised_signal, Fs_1000);
hilbert_envelope = get_hilbert_envelope(denoised_signal, Fs_1000);

shannon_envelope  = Shannon_Envelope(denoised_signal,Fs_1000); 


f = @() get_psd_envelope(denoised_signal, Fs_1000);
t.psd_envelope=timeit(f);
f = @() get_wavelet_envelope(denoised_signal, Fs_1000);
t.wavelet_envelope=timeit(f);
f = @() get_homomorphic_envelope(denoised_signal, Fs_1000);
t.homomorphic_envelope=timeit(f);
f = @() get_hilbert_envelope(denoised_signal, Fs_1000);
t.hilbert_envelope=timeit(f);

f = @() Shannon_Envelope(denoised_signal,Fs_1000); 
t.shannon_envelope=timeit(f);

% Automated Signal Quality Assessment for Heart Sound Signal by Novel Features and Evaluation in Open Public Datasets
% https://downloads.hindawi.com/journals/bmri/2021/7565398.pdf
stft_envelope=getEnvelopeFromSTFT(denoised_signal,Fs_1000);

f = @() getEnvelopeFromSTFT(denoised_signal,Fs_1000);
t.stft_envelope=timeit(f);

%% Envelope based features
% Variance-----------------------------------------------------------------
% Automated Signal Quality Assessment for Heart Sound Signal by Novel Features and Evaluation in Open Public Datasets
% https://downloads.hindawi.com/journals/bmri/2021/7565398.pdf
%SQI.psd_envelope_var = get_variance_score(psd_envelope);
%SQI.psd_envelope_original_var = get_variance_score(psd_original_envelope);
%SQI.wavelet_envelope_var = get_variance_score(wavelet_envelope);
%SQI.homomorphic_envelope_var = get_variance_score(homomorphic_envelope);
%SQI.hilbert_envelope_var = get_variance_score(hilbert_envelope);
%SQI.stft_envelope_var = get_variance_score(stft_envelope);
%SQI.shannon_envelope_var = get_variance_score(shannon_envelope);

f = @() get_variance_score(psd_envelope);
t.psd_envelope_var=timeit(f);
f = @() get_variance_score(psd_original_envelope);
t.psd_original_envelope_var=timeit(f);
f = @() get_variance_score(wavelet_envelope);
t.wavelet_envelope_var=timeit(f);
f = @() get_variance_score(homomorphic_envelope);
t.homomorphic_envelope_var=timeit(f);
f = @() get_variance_score(hilbert_envelope);
t.hilbert_envelope_var=timeit(f);
f = @() get_variance_score(stft_envelope);
t.stft_envelope_var=timeit(f);
f = @() get_variance_score(shannon_envelope);
t.shannon_envelope_var=timeit(f);

% Sample entropy-----------------------------------------------------------
% Automated Signal Quality Assessment for Heart Sound Signal by Novel Features and Evaluation in Open Public Datasets
% https://downloads.hindawi.com/journals/bmri/2021/7565398.pdf
% M=2, r=0.2
psd_envelope_down=resample(psd_envelope,30,Fs_1000);
psd_original_envelope_down=resample(psd_original_envelope,30,Fs_1000);
wavelet_envelope_down=resample(wavelet_envelope,30,Fs_1000);
homomorphic_envelope_down=resample(homomorphic_envelope,30,Fs_1000);
hilbert_envelope_down=resample(hilbert_envelope,30,Fs_1000);
stft_envelope_down=resample(stft_envelope,30,Fs_1000);
shannon_envelope_down=resample(shannon_envelope,30,Fs_1000);
%SQI.psd_envelope_se = get_entropy(psd_envelope_down, 2, 0.2);
%SQI.psd_original_envelope_se = get_entropy(psd_original_envelope_down, 2, 0.2);
%SQI.wavelet_envelope_se = get_entropy(wavelet_envelope_down, 2, 0.2);
%SQI.homomorphic_envelope_se = get_entropy(homomorphic_envelope_down, 2, 0.2);
%SQI.hilbert_envelope_se = get_entropy(hilbert_envelope_down, 2, 0.2);
%SQI.stft_envelope_se = get_entropy(stft_envelope_down, 2, 0.2);
%SQI.shannon_envelope_se = get_entropy(shannon_envelope_down, 2, 0.2);

f = @() resample(psd_envelope,30,Fs_1000);
t.psd_envelope_down=timeit(f);
f = @() resample(psd_original_envelope,30,Fs_1000);
t.psd_original_envelope_down=timeit(f);
f = @() resample(wavelet_envelope,30,Fs_1000);
t.wavelet_envelope_down=timeit(f);
f = @() resample(homomorphic_envelope,30,Fs_1000);
t.homomorphic_envelope_down=timeit(f);
f = @() resample(hilbert_envelope,30,Fs_1000);
t.hilbert_envelope_down=timeit(f);
f = @() resample(stft_envelope,30,Fs_1000);
t.stft_envelope_down=timeit(f);
f = @() resample(shannon_envelope,30,Fs_1000);
t.shannon_envelope_down=timeit(f);

f = @() get_entropy(psd_envelope_down, 2, 0.2);
t.psd_envelope_se=timeit(f);
f = @() get_entropy(psd_original_envelope_down, 2, 0.2);
t.psd_original_envelope_se=timeit(f);
f = @() get_entropy(wavelet_envelope_down, 2, 0.2);
t.wavelet_envelope_se=timeit(f);
f = @() get_entropy(homomorphic_envelope_down, 2, 0.2);
t.homomorphic_envelope_se=timeit(f);
f = @() get_entropy(hilbert_envelope_down, 2, 0.2);
t.hilbert_envelope_se=timeit(f);
f = @() get_entropy(stft_envelope_down, 2, 0.2);
t.stft_envelope_se=timeit(f);
f = @() get_entropy(shannon_envelope_down, 2, 0.2);
t.shannon_envelope_se=timeit(f);

% Average and std of maximum coefficients of adjactive envelope cycles-----
% Automated Signal Quality Assessment for Heart Sound Signal by Novel Features and Evaluation in Open Public Datasets
% https://downloads.hindawi.com/journals/bmri/2021/7565398.pdf
%[SQI.psd_envelope_cycle_corr_avg, SQI.psd_envelope_cycle_corr_std]=...
%    getAverCoefStdCoef_modified(psd_envelope, SQI.autocorr_cycle_dur_heart);
%[SQI.psd_original_envelope_cycle_corr_avg, SQI.psd_original_envelope_cycle_corr_std]=...
%    getAverCoefStdCoef_modified(psd_original_envelope, SQI.autocorr_cycle_dur_heart);
%[SQI.wavelet_envelope_cycle_corr_avg, SQI.wavelet_envelope_cycle_corr_std]=...
%    getAverCoefStdCoef_modified(wavelet_envelope, SQI.autocorr_cycle_dur_heart);
%[SQI.homomorphic_envelope_cycle_corr_avg, SQI.homomorphic_envelope_cycle_corr_std]=...
%    getAverCoefStdCoef_modified(homomorphic_envelope, SQI.autocorr_cycle_dur_heart);
%[SQI.hilbert_envelope_cycle_corr_avg, SQI.hilbert_envelope_cycle_corr_std]=...
%    getAverCoefStdCoef_modified(hilbert_envelope, SQI.autocorr_cycle_dur_heart);
%[SQI.stft_envelope_cycle_corr_avg, SQI.stft_envelope_cycle_corr_std]=...
%    getAverCoefStdCoef_modified(stft_envelope, SQI.autocorr_cycle_dur_heart);
%[SQI.shannon_envelope_cycle_corr_avg, SQI.shannon_envelope_cycle_corr_std]=...
%    getAverCoefStdCoef_modified(shannon_envelope, SQI.autocorr_cycle_dur_heart);

f = @() getAverCoefStdCoef_modified(psd_envelope, SQI.autocorr_cycle_dur_heart);
t.psd_envelope_cycle=timeit(f);
f = @() getAverCoefStdCoef_modified(psd_original_envelope, SQI.autocorr_cycle_dur_heart);
t.psd_original_envelope_cycle=timeit(f);
f = @() getAverCoefStdCoef_modified(wavelet_envelope, SQI.autocorr_cycle_dur_heart);
t.wavelet_envelope_cycle=timeit(f);
f = @() getAverCoefStdCoef_modified(homomorphic_envelope, SQI.autocorr_cycle_dur_heart);
t.homomorphic_envelope_cycle=timeit(f);
f = @() getAverCoefStdCoef_modified(hilbert_envelope, SQI.autocorr_cycle_dur_heart);
t.hilbert_envelope_cycle=timeit(f);
f = @()  getAverCoefStdCoef_modified(stft_envelope, SQI.autocorr_cycle_dur_heart);
t.stft_envelope_cycle=timeit(f);
f = @()  getAverCoefStdCoef_modified(shannon_envelope, SQI.autocorr_cycle_dur_heart);
t.shannon_envelope_cycle=timeit(f);

% std of time varying heart rate-------------------------------------------
% Automated Signal Quality Assessment for Heart Sound Signal by Novel Features and Evaluation in Open Public Datasets
% https://downloads.hindawi.com/journals/bmri/2021/7565398.pdf
%[~,SQI.psd_envelope_hr_std]=getStdHeartRate_modified(psd_envelope,max_HR,min_HR,Fs_1000);
%[~,SQI.psd_original_envelope_hr_std]=getStdHeartRate_modified(psd_original_envelope,max_HR,min_HR,Fs_1000);
%[~,SQI.wavelet_envelope_hr_std]=getStdHeartRate_modified(wavelet_envelope,max_HR,min_HR,Fs_1000);
%[~,SQI.homomorphic_envelope_hr_std]=getStdHeartRate_modified(homomorphic_envelope,max_HR,min_HR,Fs_1000);
%[~,SQI.hilbert_envelope_hr_std]=getStdHeartRate_modified(hilbert_envelope,max_HR,min_HR,Fs_1000);
%[~,SQI.stft_envelope_hr_std]=getStdHeartRate_modified(stft_envelope,max_HR,min_HR,Fs_1000);
%[~,SQI.shannon_envelope_hr_std]=getStdHeartRate_modified(shannon_envelope,max_HR,min_HR,Fs_1000);

f = @() getStdHeartRate_modified(psd_envelope,max_HR,min_HR,Fs_1000);
t.psd_envelope_hr_std=timeit(f);
f = @() getStdHeartRate_modified(psd_original_envelope,max_HR,min_HR,Fs_1000);
t.psd_original_envelope_hr_std=timeit(f);
f = @() getStdHeartRate_modified(wavelet_envelope,max_HR,min_HR,Fs_1000);
t.wavelet_envelope_hr_std=timeit(f);
f = @() getStdHeartRate_modified(homomorphic_envelope,max_HR,min_HR,Fs_1000);
t.homomorphic_envelope_hr_std=timeit(f);
f = @() getStdHeartRate_modified(hilbert_envelope,max_HR,min_HR,Fs_1000);
t.hilbert_envelope_hr_std=timeit(f);
f = @() getStdHeartRate_modified(stft_envelope,max_HR,min_HR,Fs_1000);
t.stft_envelope_hr_std=timeit(f);
f = @() getStdHeartRate_modified(shannon_envelope,max_HR,min_HR,Fs_1000);
t.shannon_envelope_hr_std=timeit(f);
%% Segmentation derieved features
% Agreement between segmentaton methods------------------------------------
% Automated signal quality assessment of mobile phone-recorded heart sound signals
% https://www.tandfonline.com/doi/pdf/10.1080/03091902.2016.1213902?needAccess=true;
%SQI.seg_agreement1 = get_bSQI_modified(peak_pos_schmidt,peak_pos_liang,Fs_1000);
%SQI.seg_agreement2 = get_bSQI_modified(peak_pos_springer,peak_pos_liang,Fs_1000);
%SQI.seg_agreement3 = get_bSQI_modified(peak_pos_springer,peak_pos_schmidt,Fs_1000);

f = @() get_bSQI_modified(peak_pos_schmidt,peak_pos_liang,Fs_1000);
t.seg_agreement1=timeit(f);
f = @() get_bSQI_modified(peak_pos_springer,peak_pos_liang,Fs_1000);
t.seg_agreement2=timeit(f);
f = @() get_bSQI_modified(peak_pos_springer,peak_pos_schmidt,Fs_1000);
t.seg_agreement3=timeit(f);

% S1 and S2 quality factor-------------------------------------------------
% Automatic signal quality index determination of radar-recorded heart sound signals using ensemble classification
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8731709
%[SQI.seg_hsmm_qf_s1_psd,SQI.seg_hsmm_qf_s2_psd] = get_hsmm_qf(assigned_states,psd_envelope);
%[SQI.seg_hsmm_qf_s1_psd_ori,SQI.seg_hsmm_qf_s2_psd_ori] = get_hsmm_qf(assigned_states,psd_original_envelope);
%[SQI.seg_hsmm_qf_s1_wav,SQI.seg_hsmm_qf_s2_wav] = get_hsmm_qf(assigned_states,wavelet_envelope);
%[SQI.seg_hsmm_qf_s1_homo,SQI.seg_hsmm_qf_homo] = get_hsmm_qf(assigned_states,homomorphic_envelope);
%[SQI.seg_hsmm_qf_s1_hil,SQI.seg_hsmm_qf_hil] = get_hsmm_qf(assigned_states,hilbert_envelope);
%[SQI.seg_hsmm_qf_s1_stft,SQI.seg_hsmm_qf_stft] = get_hsmm_qf(assigned_states,stft_envelope);
%[SQI.seg_hsmm_qf_s1_shannon,SQI.seg_hsmm_qf_stft] = get_hsmm_qf(assigned_states,shannon_envelope);

f = @() get_hsmm_qf(assigned_states,psd_envelope);
t.psd_envelope_hsmm_qf=timeit(f);
f = @() get_hsmm_qf(assigned_states,psd_original_envelope);
t.psd_original_envelope_hsmm_qf=timeit(f);
f = @() get_hsmm_qf(assigned_states,wavelet_envelope);
t.wavelet_envelope_hsmm_qf=timeit(f);
f = @() get_hsmm_qf(assigned_states,homomorphic_envelope);
t.homomorphic_envelope_hsmm_qf=timeit(f);
f = @() get_hsmm_qf(assigned_states,hilbert_envelope);
t.hilbert_envelope_hsmm_qf=timeit(f);
f = @() get_hsmm_qf(assigned_states,stft_envelope);
t.stft_envelope_hsmm_qf=timeit(f);
f = @() get_hsmm_qf(assigned_states,shannon_envelope);
t.shannon_envelope_hsmm_qf=timeit(f);

% Percentage abnormal windows----------------------------------------------
% PCG classification using a neural network approach
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7868946
%[SQI.seg_per_abnormal_rmssd,SQI.seg_per_abnormal_zero,SQI.seg_per_abnormal_sd1] = ...
%    get_segmentation_features(assigned_states,audio_downsampled_1000);

f = @() get_segmentation_features(assigned_states,audio_downsampled_1000);
t.seg_per_abnormal=timeit(f);

% Linear dependency of rows------------------------------------------------
% Noise detection during heart sound recording using periodicity signatures
% https://www.researchgate.net/publication/51037539_Noise_detection_during_heart_sound_recording_using_periodicity_signatures
%[SQI.psd_rho1,SQI.psd_rho2,SQI.psd_rho3,SQI.psd_noise] = linear_depedency_psd(psd_full,psd_full_freq);

f = @() linear_depedency_psd(psd_full,psd_full_freq);
t.seg_depedency=timeit(f);
%% Get MFCC
% Improving the quality of point of care diagnostics with real-time machine learning in low literacy LMIC settings
% https://dl.acm.org/doi/pdf/10.1145/3209811.3209815
%13 features
%Frame the window using a hamming window and a window length of 25 ms with a sliding window of 10 ms.
coeffs_detailed = mfcc(audio_data,Fs,'WindowLength',round(0.025*Fs),...
    'OverlapLength',round(0.015*Fs),'NumCoeffs',13);

f = @() mfcc(audio_data,Fs,'WindowLength',round(0.025*Fs),...
    'OverlapLength',round(0.015*Fs),'NumCoeffs',13);
t.mfcc_coeffs_detailed=timeit(f);

% Automatic signal quality index determination of radar-recorded heart sound signals using ensemble classification
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8731709
%SQI.mfcc_coeffs = mfcc(audio_data,Fs,'WindowLength',length(audio_data),'NumCoeffs',13);

f = @() mfcc(audio_data,Fs,'WindowLength',length(audio_data),'NumCoeffs',13);
t.mfcc_coeffs=timeit(f);

%% MFCC derived features
% Improving the quality of point of care diagnostics with real-time machine learning in low literacy LMIC settings
% https://dl.acm.org/doi/pdf/10.1145/3209811.3209815
% Heart sound anomaly and quality detection using ensemble of neural networks without segmentation
% https://ieeexplore.ieee.org/document/7868817
%[SQI.mfcc_min_coeffs,SQI.mfcc_max_coeffs,SQI.mfcc_mean_coeffs,SQI.mfcc_median_coeffs,...
%    SQI.mfcc_mode_coeffs,SQI.mfcc_var_coeffs,SQI.mfcc_mean_min,SQI.mfcc_mean_max,SQI.mfcc_mean_skew,SQI.mfcc_skew_coeffs]...
%    =get_mfcc_features(coeffs_detailed);

f = @() get_mfcc_features(coeffs_detailed);
t.mfcc_features=timeit(f);

%% Iga Wavelet
% PCG classification using a neural network approach
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7868946
PCG=resample(audio_data,1000,Fs);
[PCG_Hi_Hi] = wavelet_transform_iga(PCG,1000);
[qrs_pos,PCG_resampled]=checkingZero2_modified(PCG_Hi_Hi,1000);

f = @() resample(audio_data,1000,Fs);
t.wav_iga1=timeit(f);
f = @() wavelet_transform_iga(PCG,1000);
t.wav_iga2=timeit(f);
f = @() checkingZero2_modified(PCG_Hi_Hi,1000);
t.wav_iga3=timeit(f);

%% Iga Wavelet derived features
% PCG classification using a neural network approach
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7868946
%[SQI.wav_iga_per_peak_window_original,SQI.wav_iga_per_peak_window_modified,...
%    SQI.wav_iga_per_peak_window_modified2]= acceptable_window_peaks(PCG_resampled,qrs_pos,1000);
%SQI.wav_iga_transitCounterNorm =get_zcr_1(PCG_resampled);
%SQI.wav_iga_reducedTransitCounter=get_zcr_2(PCG_resampled, qrs_pos);
%SQI.wav_iga_RMSSDall=get_rmss(PCG_resampled);

f = @() acceptable_window_peaks(PCG_resampled,qrs_pos,1000);
t.wav_iga_per_peak_window=timeit(f);
f = @() get_zcr_1(PCG_resampled);
t.wav_iga_transitCounterNorm=timeit(f);
f = @() get_zcr_2(PCG_resampled, qrs_pos);
t.wav_iga_reducedTransitCounter=timeit(f);
f = @() get_rmss(PCG_resampled);
t.wav_iga_RMSSDall=timeit(f);

end






