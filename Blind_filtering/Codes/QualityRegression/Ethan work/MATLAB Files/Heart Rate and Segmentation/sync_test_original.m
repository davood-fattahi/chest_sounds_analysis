HR_est=table;
window_length=3;
%% Liang
options.env='none';
%'none','hilbert','homomorphic','psd','wavelet','stft','shannon'
options.autocorr='none';
%'none','original','filtered'
options.init_hr='none';
% 'none','autocorr_peak','autocorr_bestpeak','autocorr_findpeaks','autocorr_freq','autocorr_cc','autocorr_svd','envelope_freq', 'envelope_findpeaks','periodicity'
options.systolic='none';
%'none','yes'
options.seg='liang';
% 'none','liang''schmidt','springer'
max_HR=220;
min_HR=70; 

[HR_test,~,Comparison]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length);

HR_est.liang_peak=HR_test.seg_numpeaks/2;
HR_est.liang_hr=HR_test.seg_hr; 

%% Periodicity
options.env='none';
%'none','hilbert','homomorphic','psd','wavelet','stft','shannon'
options.autocorr='none';
%'none','original','filtered'
options.init_hr='periodicity';
% 'none','autocorr_peak','autocorr_bestpeak','autocorr_findpeaks','autocorr_freq','autocorr_cc','autocorr_svd','envelope_freq', 'envelope_findpeaks','periodicity'
options.systolic='none';
%'none','yes'
options.seg='none';
% 'none','liang''schmidt','springer'

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.periodicity_bestpeak_hr=HR_test.initial;

%% Schmidt
options.env='homomorphic';
%'none','hilbert','homomorphic','psd','wavelet','stft','shannon'
options.autocorr='none';
%'none','original','filtered'

options.init_hr='envelope_freq';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.homo_envelope_freq_hr=HR_test.initial;

options.init_hr='envelope_findpeaks';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.homo_envelope_findpeaks_hr=HR_test.initial;
HR_est.homo_envelope_findpeaks_peak=HR_test.num_peaks/2;


options.autocorr='original';
%'none','original','filtered'
options.init_hr='autocorr_peak';
% 'none','autocorr_peak','autocorr_bestpeak','autocorr_findpeaks','autocorr_freq','autocorr_cc','autocorr_svd','envelope_freq', 'envelope_findpeaks','periodicity'

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.homo_ori_autocorr_peak_hr=HR_test.initial;

options.init_hr='autocorr_bestpeak';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.homo_ori_autocorr_bestpeak_hr=HR_test.initial;

options.init_hr='autocorr_findpeaks';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.homo_ori_autocorr_findpeaks_hr=HR_test.initial;
HR_est.homo_ori_autocorr_findpeaks_peak=HR_test.num_peaks;

options.init_hr='autocorr_freq';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.homo_ori_autocorr_freq_hr=HR_test.initial;

options.init_hr='autocorr_cc';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.homo_ori_autocorr_cc_hr=HR_test.initial;

options.init_hr='autocorr_svd';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.homo_ori_autocorr_svd_hr=HR_test.initial;



options.autocorr='filtered';
options.init_hr='autocorr_peak';
% 'none','autocorr_peak','autocorr_bestpeak','autocorr_findpeaks','autocorr_freq','autocorr_cc','autocorr_svd','envelope_freq', 'envelope_findpeaks','periodicity'

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.homo_fil_autocorr_peak_hr=HR_test.initial;

options.init_hr='autocorr_bestpeak';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.homo_fil_autocorr_bestpeak_hr=HR_test.initial;

options.init_hr='autocorr_findpeaks';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.homo_fil_autocorr_findpeaks_hr=HR_test.initial;
HR_est.homo_fil_autocorr_findpeaks_peak=HR_test.num_peaks;

options.init_hr='autocorr_freq';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.homo_fil_autocorr_freq_hr=HR_test.initial;

options.init_hr='autocorr_cc';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.homo_fil_autocorr_cc_hr=HR_test.initial;

options.init_hr='autocorr_svd';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.homo_fil_autocorr_svd_hr=HR_test.initial;















%% Hilbert
options.env='hilbert';
%'none','hilbert','homomorphic','psd','wavelet','stft','shannon'
options.autocorr='none';
%'none','original','filtered'

options.init_hr='envelope_freq';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.hil_envelope_freq_hr=HR_test.initial;

options.init_hr='envelope_findpeaks';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.hil_envelope_findpeaks_hr=HR_test.initial;
HR_est.hil_envelope_findpeaks_peak=HR_test.num_peaks/2;


options.autocorr='original';
%'none','original','filtered'
options.init_hr='autocorr_peak';
% 'none','autocorr_peak','autocorr_bestpeak','autocorr_findpeaks','autocorr_freq','autocorr_cc','autocorr_svd','envelope_freq', 'envelope_findpeaks','periodicity'

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.hil_ori_autocorr_peak_hr=HR_test.initial;

options.init_hr='autocorr_bestpeak';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.hil_ori_autocorr_bestpeak_hr=HR_test.initial;

options.init_hr='autocorr_findpeaks';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.hil_ori_autocorr_findpeaks_hr=HR_test.initial;
HR_est.hil_ori_autocorr_findpeaks_peak=HR_test.num_peaks;

options.init_hr='autocorr_freq';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.hil_ori_autocorr_freq_hr=HR_test.initial;

options.init_hr='autocorr_cc';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.hil_ori_autocorr_cc_hr=HR_test.initial;

options.init_hr='autocorr_svd';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.hil_ori_autocorr_svd_hr=HR_test.initial;



options.autocorr='filtered';
options.init_hr='autocorr_peak';
% 'none','autocorr_peak','autocorr_bestpeak','autocorr_findpeaks','autocorr_freq','autocorr_cc','autocorr_svd','envelope_freq', 'envelope_findpeaks','periodicity'

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.hil_fil_autocorr_peak_hr=HR_test.initial;

options.init_hr='autocorr_bestpeak';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.hil_fil_autocorr_bestpeak_hr=HR_test.initial;

options.init_hr='autocorr_findpeaks';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.hil_fil_autocorr_findpeaks_hr=HR_test.initial;
HR_est.hil_fil_autocorr_findpeaks_peak=HR_test.num_peaks;

options.init_hr='autocorr_freq';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.hil_fil_autocorr_freq_hr=HR_test.initial;

options.init_hr='autocorr_cc';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.hil_fil_autocorr_cc_hr=HR_test.initial;

options.init_hr='autocorr_svd';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.hil_fil_autocorr_svd_hr=HR_test.initial;













%% Schmidt
options.env='psd';
%'none','hilbert','homomorphic','psd','wavelet','stft','shannon'
options.autocorr='none';
%'none','original','filtered'

options.init_hr='envelope_freq';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.psd_envelope_freq_hr=HR_test.initial;

options.init_hr='envelope_findpeaks';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.psd_envelope_findpeaks_hr=HR_test.initial;
HR_est.psd_envelope_findpeaks_peak=HR_test.num_peaks/2;


options.autocorr='original';
%'none','original','filtered'
options.init_hr='autocorr_peak';
% 'none','autocorr_peak','autocorr_bestpeak','autocorr_findpeaks','autocorr_freq','autocorr_cc','autocorr_svd','envelope_freq', 'envelope_findpeaks','periodicity'

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.psd_ori_autocorr_peak_hr=HR_test.initial;

options.init_hr='autocorr_bestpeak';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.psd_ori_autocorr_bestpeak_hr=HR_test.initial;

options.init_hr='autocorr_findpeaks';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.psd_ori_autocorr_findpeaks_hr=HR_test.initial;
HR_est.psd_ori_autocorr_findpeaks_peak=HR_test.num_peaks;

options.init_hr='autocorr_freq';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.psd_ori_autocorr_freq_hr=HR_test.initial;

options.init_hr='autocorr_cc';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.psd_ori_autocorr_cc_hr=HR_test.initial;

options.init_hr='autocorr_svd';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.psd_ori_autocorr_svd_hr=HR_test.initial;


options.autocorr='filtered';
options.init_hr='autocorr_peak';
% 'none','autocorr_peak','autocorr_bestpeak','autocorr_findpeaks','autocorr_freq','autocorr_cc','autocorr_svd','envelope_freq', 'envelope_findpeaks','periodicity'

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.psd_fil_autocorr_peak_hr=HR_test.initial;

options.init_hr='autocorr_bestpeak';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.psd_fil_autocorr_bestpeak_hr=HR_test.initial;

options.init_hr='autocorr_findpeaks';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.psd_fil_autocorr_findpeaks_hr=HR_test.initial;
HR_est.psd_fil_autocorr_findpeaks_peak=HR_test.num_peaks;

options.init_hr='autocorr_freq';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.psd_fil_autocorr_freq_hr=HR_test.initial;

options.init_hr='autocorr_cc';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.psd_fil_autocorr_cc_hr=HR_test.initial;

options.init_hr='autocorr_svd';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.psd_fil_autocorr_svd_hr=HR_test.initial;











%% Schmidt
options.env='wavelet';
%'none','hilbert','homomorphic','psd','wavelet','stft','shannon'
options.autocorr='none';
%'none','original','filtered'

options.init_hr='envelope_freq';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.wav_envelope_freq_hr=HR_test.initial;

options.init_hr='envelope_findpeaks';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.wav_envelope_findpeaks_hr=HR_test.initial;
HR_est.wav_envelope_findpeaks_peak=HR_test.num_peaks/2;


options.autocorr='original';
%'none','original','filtered'
options.init_hr='autocorr_peak';
% 'none','autocorr_peak','autocorr_bestpeak','autocorr_findpeaks','autocorr_freq','autocorr_cc','autocorr_svd','envelope_freq', 'envelope_findpeaks','periodicity'

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.wav_ori_autocorr_peak_hr=HR_test.initial;

options.init_hr='autocorr_bestpeak';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.wav_ori_autocorr_bestpeak_hr=HR_test.initial;

options.init_hr='autocorr_findpeaks';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.wav_ori_autocorr_findpeaks_hr=HR_test.initial;
HR_est.wav_ori_autocorr_findpeaks_peak=HR_test.num_peaks;

options.init_hr='autocorr_freq';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.wav_ori_autocorr_freq_hr=HR_test.initial;

options.init_hr='autocorr_cc';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.wav_ori_autocorr_cc_hr=HR_test.initial;

options.init_hr='autocorr_svd';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.wav_ori_autocorr_svd_hr=HR_test.initial;


options.autocorr='filtered';
options.init_hr='autocorr_peak';
% 'none','autocorr_peak','autocorr_bestpeak','autocorr_findpeaks','autocorr_freq','autocorr_cc','autocorr_svd','envelope_freq', 'envelope_findpeaks','periodicity'

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.wav_fil_autocorr_peak_hr=HR_test.initial;

options.init_hr='autocorr_bestpeak';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.wav_fil_autocorr_bestpeak_hr=HR_test.initial;

options.init_hr='autocorr_findpeaks';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.wav_fil_autocorr_findpeaks_hr=HR_test.initial;
HR_est.wav_fil_autocorr_findpeaks_peak=HR_test.num_peaks;

options.init_hr='autocorr_freq';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.wav_fil_autocorr_freq_hr=HR_test.initial;

options.init_hr='autocorr_cc';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.wav_fil_autocorr_cc_hr=HR_test.initial;

options.init_hr='autocorr_svd';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.wav_fil_autocorr_svd_hr=HR_test.initial;











%% Schmidt
options.env='stft';
%'none','hilbert','homomorphic','psd','wavelet','stft','shannon'
options.autocorr='none';
%'none','original','filtered'

options.init_hr='envelope_freq';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.stft_envelope_freq_hr=HR_test.initial;

options.init_hr='envelope_findpeaks';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.stft_envelope_findpeaks_hr=HR_test.initial;
HR_est.stft_envelope_findpeaks_peak=HR_test.num_peaks/2;


options.autocorr='original';
%'none','original','filtered'
options.init_hr='autocorr_peak';
% 'none','autocorr_peak','autocorr_bestpeak','autocorr_findpeaks','autocorr_freq','autocorr_cc','autocorr_svd','envelope_freq', 'envelope_findpeaks','periodicity'

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.stft_ori_autocorr_peak_hr=HR_test.initial;

options.init_hr='autocorr_bestpeak';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.stft_ori_autocorr_bestpeak_hr=HR_test.initial;

options.init_hr='autocorr_findpeaks';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.stft_ori_autocorr_findpeaks_hr=HR_test.initial;
HR_est.stft_ori_autocorr_findpeaks_peak=HR_test.num_peaks;

options.init_hr='autocorr_freq';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.stft_ori_autocorr_freq_hr=HR_test.initial;

options.init_hr='autocorr_cc';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.stft_ori_autocorr_cc_hr=HR_test.initial;

options.init_hr='autocorr_svd';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.stft_ori_autocorr_svd_hr=HR_test.initial;


options.autocorr='filtered';
options.init_hr='autocorr_peak';
% 'none','autocorr_peak','autocorr_bestpeak','autocorr_findpeaks','autocorr_freq','autocorr_cc','autocorr_svd','envelope_freq', 'envelope_findpeaks','periodicity'

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.stft_fil_autocorr_peak_hr=HR_test.initial;

options.init_hr='autocorr_bestpeak';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.stft_fil_autocorr_bestpeak_hr=HR_test.initial;

options.init_hr='autocorr_findpeaks';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.stft_fil_autocorr_findpeaks_hr=HR_test.initial;
HR_est.stft_fil_autocorr_findpeaks_peak=HR_test.num_peaks;

options.init_hr='autocorr_freq';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.stft_fil_autocorr_freq_hr=HR_test.initial;

options.init_hr='autocorr_cc';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.stft_fil_autocorr_cc_hr=HR_test.initial;

options.init_hr='autocorr_svd';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.stft_fil_autocorr_svd_hr=HR_test.initial;










%% Schmidt
options.env='shannon';
%'none','hilbert','homomorphic','psd','wavelet','stft','shannon'
options.autocorr='none';
%'none','original','filtered'

options.init_hr='envelope_freq';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.sha_envelope_freq_hr=HR_test.initial;

options.init_hr='envelope_findpeaks';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.sha_envelope_findpeaks_hr=HR_test.initial;
HR_est.sha_envelope_findpeaks_peak=HR_test.num_peaks/2;


options.autocorr='original';
%'none','original','filtered'
options.init_hr='autocorr_peak';
% 'none','autocorr_peak','autocorr_bestpeak','autocorr_findpeaks','autocorr_freq','autocorr_cc','autocorr_svd','envelope_freq', 'envelope_findpeaks','periodicity'

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.sha_ori_autocorr_peak_hr=HR_test.initial;

options.init_hr='autocorr_bestpeak';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.sha_ori_autocorr_bestpeak_hr=HR_test.initial;

options.init_hr='autocorr_findpeaks';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.sha_ori_autocorr_findpeaks_hr=HR_test.initial;
HR_est.sha_ori_autocorr_findpeaks_peak=HR_test.num_peaks;

options.init_hr='autocorr_freq';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.sha_ori_autocorr_freq_hr=HR_test.initial;

options.init_hr='autocorr_cc';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.sha_ori_autocorr_cc_hr=HR_test.initial;

options.init_hr='autocorr_svd';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.sha_ori_autocorr_svd_hr=HR_test.initial;


options.autocorr='filtered';
options.init_hr='autocorr_peak';
% 'none','autocorr_peak','autocorr_bestpeak','autocorr_findpeaks','autocorr_freq','autocorr_cc','autocorr_svd','envelope_freq', 'envelope_findpeaks','periodicity'

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.sha_fil_autocorr_peak_hr=HR_test.initial;

options.init_hr='autocorr_bestpeak';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.sha_fil_autocorr_bestpeak_hr=HR_test.initial;

options.init_hr='autocorr_findpeaks';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.sha_fil_autocorr_findpeaks_hr=HR_test.initial;
HR_est.sha_fil_autocorr_findpeaks_peak=HR_test.num_peaks;

options.init_hr='autocorr_freq';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.sha_fil_autocorr_freq_hr=HR_test.initial;

options.init_hr='autocorr_cc';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.sha_fil_autocorr_cc_hr=HR_test.initial;

options.init_hr='autocorr_svd';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length); 
HR_est.sha_fil_autocorr_svd_hr=HR_test.initial;





options.env='shannon';
%'none','hilbert','homomorphic','psd','wavelet','stft','shannon'
options.autocorr='filtered';
%'none','original','filtered'
options.init_hr='autocorr_bestpeak';
% 'none','autocorr_peak','autocorr_bestpeak','autocorr_findpeaks','autocorr_freq','autocorr_cc','autocorr_svd','envelope_freq', 'envelope_findpeaks','periodicity'
options.systolic='yes';
%'none','yes'
options.seg='springer';
% 'none','liang''schmidt','springer'

load('Springer_B_matrix.mat', 'Springer_B_matrix');
load('Springer_pi_vector.mat', 'Springer_pi_vector');
load('Springer_total_obs_distribution.mat', 'Springer_total_obs_distribution');

options.seg_fs=50;
options.seg_pi_vector=Springer_pi_vector;
options.seg_b_matrix=Springer_B_matrix;
options.seg_total_obs_dist=Springer_total_obs_distribution;


[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length);

HR_est.sha_fil_autocorr_bestpeak_springer_hr=HR_test.seg_hr;
HR_est.sha_fil_autocorr_bestpeak_springer_peak=HR_test.seg_numpeaks;

options.init_hr='autocorr_peak';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length);

HR_est.sha_fil_autocorr_peak_springer_hr=HR_test.seg_hr;
HR_est.sha_fil_autocorr_peak_springer_peak=HR_test.seg_numpeaks;



%'none','hilbert','homomorphic','psd','wavelet','stft','shannon'
options.env='hilbert';
options.init_hr='autocorr_peak';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length);

HR_est.hil_fil_autocorr_bestpeak_springer_hr=HR_test.seg_hr;
HR_est.hil_fil_autocorr_bestpeak_springer_peak=HR_test.seg_numpeaks;

options.init_hr='autocorr_peak';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length);

HR_est.hil_fil_autocorr_peak_springer_hr=HR_test.seg_hr;
HR_est.hil_fil_autocorr_peak_springer_peak=HR_test.seg_numpeaks;



%'none','hilbert','homomorphic','psd','wavelet','stft','shannon'
options.env='stft';
options.init_hr='autocorr_peak';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length);

HR_est.stft_fil_autocorr_bestpeak_springer_hr=HR_test.seg_hr;
HR_est.stft_fil_autocorr_bestpeak_springer_peak=HR_test.seg_numpeaks;

options.init_hr='autocorr_peak';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length);

HR_est.stft_fil_autocorr_peak_springer_hr=HR_test.seg_hr;
HR_est.stft_fil_autocorr_peak_springer_peak=HR_test.seg_numpeaks;


%'none','hilbert','homomorphic','psd','wavelet','stft','shannon'
options.env='homomorphic';
options.init_hr='autocorr_peak';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length);

HR_est.homo_fil_autocorr_bestpeak_springer_hr=HR_test.seg_hr;
HR_est.homo_fil_autocorr_bestpeak_springer_peak=HR_test.seg_numpeaks;

options.init_hr='autocorr_peak';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length);

HR_est.homo_fil_autocorr_peak_springer_hr=HR_test.seg_hr;
HR_est.homo_fil_autocorr_peak_springer_peak=HR_test.seg_numpeaks;













options.env='shannon';
%'none','hilbert','homomorphic','psd','wavelet','stft','shannon'
options.autocorr='filtered';
%'none','original','filtered'
options.init_hr='autocorr_bestpeak';
% 'none','autocorr_peak','autocorr_bestpeak','autocorr_findpeaks','autocorr_freq','autocorr_cc','autocorr_svd','envelope_freq', 'envelope_findpeaks','periodicity'
options.systolic='yes';
%'none','yes'
options.seg='schmidt';
% 'none','liang''schmidt','springer'

load('hmm.mat','b_matrix');
load('hmm.mat', 'pi_vector');

options.seg_fs=50;
options.seg_pi_vector=pi_vector;
options.seg_b_matrix=b_matrix;
options.seg_total_obs_dist=Springer_total_obs_distribution;


[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length);

HR_est.sha_fil_autocorr_bestpeak_schmidt_hr=HR_test.seg_hr;
HR_est.sha_fil_autocorr_bestpeak_schmidt_peak=HR_test.seg_numpeaks;

options.init_hr='autocorr_peak';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length);

HR_est.sha_fil_autocorr_peak_schmidt_hr=HR_test.seg_hr;
HR_est.sha_fil_autocorr_peak_schmidt_peak=HR_test.seg_numpeaks;



%'none','hilbert','homomorphic','psd','wavelet','stft','shannon'
options.env='hilbert';
options.init_hr='autocorr_peak';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length);

HR_est.hil_fil_autocorr_bestpeak_schmidt_hr=HR_test.seg_hr;
HR_est.hil_fil_autocorr_bestpeak_schmidt_peak=HR_test.seg_numpeaks;

options.init_hr='autocorr_peak';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length);

HR_est.hil_fil_autocorr_peak_schmidt_hr=HR_test.seg_hr;
HR_est.hil_fil_autocorr_peak_schmidt_peak=HR_test.seg_numpeaks;



%'none','hilbert','homomorphic','psd','wavelet','stft','shannon'
options.env='stft';
options.init_hr='autocorr_peak';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length);

HR_est.stft_fil_autocorr_bestpeak_schmidt_hr=HR_test.seg_hr;
HR_est.stft_fil_autocorr_bestpeak_schmidt_peak=HR_test.seg_numpeaks;

options.init_hr='autocorr_peak';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length);

HR_est.stft_fil_autocorr_peak_schmidt_hr=HR_test.seg_hr;
HR_est.stft_fil_autocorr_peak_schmidt_peak=HR_test.seg_numpeaks;


%'none','hilbert','homomorphic','psd','wavelet','stft','shannon'
options.env='homomorphic';
options.init_hr='autocorr_peak';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length);

HR_est.homo_fil_autocorr_bestpeak_schmidt_hr=HR_test.seg_hr; 
HR_est.homo_fil_autocorr_bestpeak_schmidt_peak=HR_test.seg_numpeaks;

options.init_hr='autocorr_peak';

[HR_test,~,~]=get_hr_test_sync(Sync_Data,File_Info,options,max_HR,min_HR,window_length);

HR_est.homo_fil_autocorr_peak_schmidt_hr=HR_test.seg_hr;
HR_est.homo_fil_autocorr_peak_schmidt_peak=HR_test.seg_numpeaks;

















 