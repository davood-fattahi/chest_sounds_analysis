clear
% close all
clc

Heart_Comparison=readtable('Sync Vital Signs.xlsx');
files=unique(Heart_Comparison.File); 

for i=1:length(files)
    file=files{i};
    locs=strcmp(Heart_Comparison.File,file); 
    BR_sec(i,1)=mean(Heart_Comparison.BR_sec(locs));
    HR_sec(i,1)=mean(Heart_Comparison.HR_sec(locs)) ;
end

%% Lung Parameters
BR_est=table; 
min_BR=15;
max_BR=100; 

%% Lung Power
options.env='psd1';
options.autocorr='filtered';
options.init_hr= 'envelope_findpeaks_br'; 
options.seg='none';

DA='..\..\ChestRecords\CurrentDataset\';
dataSheet=readcell([DA '00dataSheet.xlsx']);
files=dataSheet([false; cell2mat(dataSheet(2:end,6))==1],1);
files = cellfun(@(x) x(regexp(x,'sync_')+5:regexp(x,'.wav')-1), files,'UniformOutput',false);
file_ending=''; 
[BR_test]=get_br_test_sync_correction_10(files,options,max_BR,min_BR, DA, file_ending);
lung_raw=BR_test.initial_2;
% lung_raw(isnan(lung_raw))=0; 

DA='..\FilteringMethods\FrequencyFiltering\Results\FreqFixLung\';
file_ending='_FreqFilt_Lung'; 
[BR_test]=get_br_test_sync_correction_10(files,options,max_BR,min_BR, DA, file_ending);
lung_freqfilt=BR_test.initial_2;


DA='..\FilteringMethods\swt-pca\TempResults\SwtPcaLung\';
file_ending='_SwtPca_Lung'; 
[BR_test]=get_br_test_sync_correction_10(files,options,max_BR,min_BR, DA, file_ending);
lung_pca_swt=BR_test.initial_2;


DA='..\FilteringMethods\cwt-pca\TempResults\CwtPcaLung\';
file_ending='_CwtPca_Lung'; 
[BR_test]=get_br_test_sync_correction_10(files,options,max_BR,min_BR, DA, file_ending);
lung_pca_cwt=BR_test.initial_2;


DA='..\FilteringMethods\cwt-pica\TempResults\CwtPicaLung\';
file_ending='_CwtPica_Lung'; 
[BR_test]=get_br_test_sync_correction_10(files,options,max_BR,min_BR, DA, file_ending);
lung_pica_cwt=BR_test.initial_2;


DA='..\FilteringMethods\cwt-sobi\TempResults\CwtSobiLung\';
file_ending='_CwtSobi_Lung'; 
[BR_test]=get_br_test_sync_correction_10(files,options,max_BR,min_BR, DA, file_ending);
lung_sobi_cwt=BR_test.initial_2;


DA='..\FilteringMethods\otherMethods\TempResults\emdLung\';
file_ending='_emd_Lung'; 
[BR_test]=get_br_test_sync_correction_10(files,options,max_BR,min_BR, DA, file_ending);
lung_emd=BR_test.initial_2;


DA='..\FilteringMethods\otherMethods\TempResults\eemdLung\';
file_ending='_eemd_Lung'; 
[BR_test]=get_br_test_sync_correction_10(files,options,max_BR,min_BR, DA, file_ending);
lung_eemd=BR_test.initial_2;


DA='..\FilteringMethods\otherMethods\TempResults\ceemdLung\';
file_ending='_ceemd_Lung'; 
[BR_test]=get_br_test_sync_correction_10(files,options,max_BR,min_BR, DA, file_ending);
lung_ceemd=BR_test.initial_2;

DA='..\FilteringMethods\otherMethods\TempResults\nmfcklLung\';
file_ending='_nmfckl_Lung'; 
[BR_test]=get_br_test_sync_correction_10(files,options,max_BR,min_BR, DA, file_ending);
lung_nmfc_kl=BR_test.initial_2;


DA='..\FilteringMethods\otherMethods\TempResults\nmfcGenFiltLung\';
file_ending='_nmfcGenFilt_Lung'; 
[BR_test]=get_br_test_sync_correction_10(files,options,max_BR,min_BR, DA, file_ending);
lung_nmfc_gen_filt=BR_test.initial_2;




DA='..\FilteringMethods\otherMethods\TempResults\nmfcSparseLung\';
file_ending='_nmfcSparse_Lung'; 
[BR_test]=get_br_test_sync_correction_10(files,options,max_BR,min_BR, DA, file_ending);
lung_nmfc_sparse=BR_test.initial_2;


DA='..\FilteringMethods\otherMethods\TempResults\ssaLung\';
file_ending='_ssa_Lung'; 
[BR_test]=get_br_test_sync_correction_10(files,options,max_BR,min_BR, DA, file_ending);
lung_ssa=BR_test.initial_2;
lung_ssa(lung_ssa==0)=nan;



lung_results=[lung_raw lung_freqfilt lung_pca_swt lung_pca_cwt lung_pica_cwt lung_sobi_cwt lung_ssa ...
    lung_emd lung_eemd lung_ceemd lung_nmfc_kl lung_nmfc_gen_filt  ...
    lung_nmfc_sparse]; 
lung_mae=mean(abs(lung_results-BR_sec),'omitnan');
lung_rmse=sqrt(mean((lung_results-BR_sec).^2,'omitnan'));

%% Heart Parameters
HR_est=table;
max_HR=180;
min_HR=130; 

% Heart Schmidt
options.env='hilbert';
options.autocorr='filtered';
options.init_hr='autocorr_findpeaks';
options.systolic='none';
options.seg='none';


DA='..\..\ChestRecords\CurrentDataset\';
dataSheet=readcell([DA '00dataSheet.xlsx']);
files=dataSheet([false; cell2mat(dataSheet(2:end,6))==1],1);
files = cellfun(@(x) x(regexp(x,'sync_')+5:regexp(x,'.wav')-1), files,'UniformOutput',false);
file_ending=''; 
[HR_test]=get_hr_test_sync_correction_10(files,options,max_HR,min_HR, DA, file_ending);
heart_raw=HR_test.initial;


DA='..\FilteringMethods\FrequencyFiltering\Results\FreqFixHeart\';
file_ending='_FreqFilt_Heart'; 
[HR_test]=get_hr_test_sync_correction_10(files,options,max_HR,min_HR, DA, file_ending);
heart_freqfilt=HR_test.initial;


DA='..\FilteringMethods\swt-pca\TempResults\SwtPcaHeart\';
file_ending='_SwtPca_Heart'; 
[HR_test]=get_hr_test_sync_correction_10(files,options,max_HR,min_HR, DA, file_ending);
heart_pca_swt=HR_test.initial;

DA='..\FilteringMethods\cwt-pca\TempResults\CwtPcaHeart\';
file_ending='_CwtPca_Heart'; 
[HR_test]=get_hr_test_sync_correction_10(files,options,max_HR,min_HR, DA, file_ending);
heart_pca_cwt=HR_test.initial;


DA='..\FilteringMethods\cwt-pica\TempResults\CwtPicaHeart\';
file_ending='_CwtPica_Heart'; 
[HR_test]=get_hr_test_sync_correction_10(files,options,max_HR,min_HR, DA, file_ending);
heart_pica_cwt=HR_test.initial;


DA='..\FilteringMethods\cwt-sobi\TempResults\CwtSobiHeart\';
file_ending='_CwtSobi_Heart'; 
[HR_test]=get_hr_test_sync_correction_10(files,options,max_HR,min_HR, DA, file_ending);
heart_sobi_cwt=HR_test.initial;



DA='..\FilteringMethods\otherMethods\TempResults\emdHeart\';
file_ending='_emd_Heart'; 
[HR_test]=get_hr_test_sync_correction_10(files,options,max_HR,min_HR, DA, file_ending);
heart_emd=HR_test.initial;


DA='..\FilteringMethods\otherMethods\TempResults\eemdHeart\';
file_ending='_eemd_Heart'; 
[HR_test]=get_hr_test_sync_correction_10(files,options,max_HR,min_HR, DA, file_ending);
heart_eemd=HR_test.initial;


DA='..\FilteringMethods\otherMethods\TempResults\ceemdHeart\';
file_ending='_ceemd_Heart'; 
[HR_test]=get_hr_test_sync_correction_10(files,options,max_HR,min_HR, DA, file_ending);
heart_ceemd=HR_test.initial;

DA='..\FilteringMethods\otherMethods\TempResults\nmfcklHeart\';
file_ending='_nmfckl_Heart'; 
[HR_test]=get_hr_test_sync_correction_10(files,options,max_HR,min_HR, DA, file_ending);
heart_nmfc_kl=HR_test.initial;


DA='..\FilteringMethods\otherMethods\TempResults\nmfcGenFiltHeart\';
file_ending='_nmfcGenFilt_Heart'; 
[HR_test]=get_hr_test_sync_correction_10(files,options,max_HR,min_HR, DA, file_ending);
heart_nmfc_gen_filt=HR_test.initial;




DA='..\FilteringMethods\otherMethods\TempResults\nmfcSparseHeart\';
file_ending='_nmfcSparse_Heart'; 
[HR_test]=get_hr_test_sync_correction_10(files,options,max_HR,min_HR, DA, file_ending);
heart_nmfc_sparse=HR_test.initial;


DA='..\FilteringMethods\otherMethods\TempResults\ssaHeart\';
file_ending='_ssa_Heart'; 
[HR_test]=get_hr_test_sync_correction_10(files,options,max_HR,min_HR, DA, file_ending);
heart_ssa=HR_test.initial;
heart_ssa(heart_ssa==0)=nan;



heart_results=[heart_raw heart_freqfilt heart_pca_swt heart_pca_cwt heart_pica_cwt heart_sobi_cwt heart_ssa ...
    heart_emd heart_eemd heart_ceemd heart_nmfc_kl heart_nmfc_gen_filt  ...
    heart_nmfc_sparse];

heart_mae=mean(abs(heart_results-HR_sec),'omitnan');
heart_rmse=sqrt(mean((heart_results-HR_sec).^2,'omitnan'));

%%

%% Heart Springer
options.env='hilbert';
options.autocorr='filtered';
options.init_hr='autocorr_findpeaks';
options.systolic='yes';
options.seg='springer';

load('Springer_B_matrix.mat', 'Springer_B_matrix');
load('Springer_pi_vector.mat', 'Springer_pi_vector');
load('Springer_total_obs_distribution.mat', 'Springer_total_obs_distribution');

options.seg_fs=50;
options.seg_pi_vector=Springer_pi_vector;
options.seg_b_matrix=Springer_B_matrix;
options.seg_total_obs_dist=Springer_total_obs_distribution;


DA='..\..\ChestRecords\CurrentDataset\';
dataSheet=readcell([DA '00dataSheet.xlsx']);
files=dataSheet([false; cell2mat(dataSheet(2:end,6))==1],1);
files = cellfun(@(x) x(regexp(x,'sync_')+5:regexp(x,'.wav')-1), files,'UniformOutput',false);
file_ending=''; 
[HR_test]=get_hr_test_sync_correction_10(files,options,max_HR,min_HR, DA, file_ending);
heart_raw=HR_test.seg_hr;


DA='..\FilteringMethods\FrequencyFiltering\Results\FreqFixHeart\';
file_ending='_FreqFilt_Heart'; 
[HR_test]=get_hr_test_sync_correction_10(files,options,max_HR,min_HR, DA, file_ending);
heart_freqfilt=HR_test.seg_hr;


DA='..\FilteringMethods\swt-pca\TempResults\SwtPcaHeart\';
file_ending='_SwtPca_Heart'; 
[HR_test]=get_hr_test_sync_correction_10(files,options,max_HR,min_HR, DA, file_ending);
heart_pca_swt=HR_test.seg_hr;


DA='..\FilteringMethods\cwt-pca\TempResults\CwtPcaHeart\';
file_ending='_CwtPca_Heart'; 
[HR_test]=get_hr_test_sync_correction_10(files,options,max_HR,min_HR, DA, file_ending);
heart_pca_cwt=HR_test.seg_hr;


DA='..\FilteringMethods\cwt-pica\TempResults\CwtPicaHeart\';
file_ending='_CwtPica_Heart'; 
[HR_test]=get_hr_test_sync_correction_10(files,options,max_HR,min_HR, DA, file_ending);
heart_pica_cwt=HR_test.seg_hr;


DA='..\FilteringMethods\cwt-sobi\TempResults\CwtSobiHeart\';
file_ending='_CwtSobi_Heart'; 
[HR_test]=get_hr_test_sync_correction_10(files,options,max_HR,min_HR, DA, file_ending);
heart_sobi_cwt=HR_test.seg_hr;




DA='..\FilteringMethods\otherMethods\TempResults\emdHeart\';
file_ending='_emd_Heart'; 
[HR_test]=get_hr_test_sync_correction_10(files,options,max_HR,min_HR, DA, file_ending);
heart_emd=HR_test.seg_hr;


DA='..\FilteringMethods\otherMethods\TempResults\eemdHeart\';
file_ending='_eemd_Heart'; 
[HR_test]=get_hr_test_sync_correction_10(files,options,max_HR,min_HR, DA, file_ending);
heart_eemd=HR_test.seg_hr;


DA='..\FilteringMethods\otherMethods\TempResults\ceemdHeart\';
file_ending='_ceemd_Heart'; 
[HR_test]=get_hr_test_sync_correction_10(files,options,max_HR,min_HR, DA, file_ending);
heart_ceemd=HR_test.seg_hr;

DA='..\FilteringMethods\otherMethods\TempResults\nmfcklHeart\';
file_ending='_nmfckl_Heart'; 
[HR_test]=get_hr_test_sync_correction_10(files,options,max_HR,min_HR, DA, file_ending);
heart_nmfc_kl=HR_test.seg_hr;


DA='..\FilteringMethods\otherMethods\TempResults\nmfcGenFiltHeart\';
file_ending='_nmfcGenFilt_Heart'; 
[HR_test]=get_hr_test_sync_correction_10(files,options,max_HR,min_HR, DA, file_ending);
heart_nmfc_gen_filt=HR_test.seg_hr;




DA='..\FilteringMethods\otherMethods\TempResults\nmfcSparseHeart\';
file_ending='_nmfcSparse_Heart'; 
[HR_test]=get_hr_test_sync_correction_10(files,options,max_HR,min_HR, DA, file_ending);
heart_nmfc_sparse=HR_test.seg_hr;


DA='..\FilteringMethods\otherMethods\TempResults\ssaHeart\';
file_ending='_ssa_Heart'; 
[HR_test]=get_hr_test_sync_correction_10(files,options,max_HR,min_HR, DA, file_ending);
heart_ssa=HR_test.seg_hr;
heart_ssa(heart_ssa==0)=nan;


heart_results2=[heart_raw heart_freqfilt heart_pca_swt heart_pca_cwt heart_pica_cwt heart_sobi_cwt heart_ssa ...
    heart_emd heart_eemd heart_ceemd heart_nmfc_kl heart_nmfc_gen_filt  ...
    heart_nmfc_sparse];

heart_mae2=mean(abs(heart_results2-HR_sec),'omitnan');
heart_rmse2=sqrt(mean((heart_results2-HR_sec).^2,'omitnan'));


%% save results
delete 'VitalSignsAnalysis.xlsx'

writecell([{'hr_NICU' 'hr_raw' 'hr_freqfilt' 'hr_swtpca' 'hr_cwtpca' 'hr_cwtpica' 'hr_cwtsobi' 'hr_ssa' 'hr_emd' ...
    'hr_eemd' 'hr_ceemd' 'hr_nmfcKL' 'hr_nmfcL2' 'hr_nmfcSparse'} ...
    ; num2cell([HR_sec heart_results])],'VitalSignsAnalysis.xlsx','Sheet', 'Heart_Schmidt');

writecell([{'hr_NICU' 'hr_raw' 'hr_freqfilt' 'hr_swtpca' 'hr_cwtpca' 'hr_cwtpica' 'hr_cwtsobi' 'hr_ssa' 'hr_emd' ...
    'hr_eemd' 'hr_ceemd' 'hr_nmfcKL' 'hr_nmfcL2' 'hr_nmfcSparse'} ...
    ; num2cell([HR_sec heart_results2])],'VitalSignsAnalysis.xlsx','Sheet', 'Heart_Springer');

writecell([{'br_NICU' 'br_raw' 'br_freqfilt' 'br_swtpca' 'br_cwtpca' 'br_cwtpica' 'br_cwtsobi' 'br_ssa' 'br_emd' ...
    'br_eemd' 'br_ceemd' 'br_nmfcKL' 'br_nmfcL2'   'br_nmfcSparse'} ...
    ; num2cell([BR_sec lung_results])],'VitalSignsAnalysis.xlsx','Sheet', 'Lung');

writecell([{'Heart_Springer' 'hr_raw' 'hr_freqfilt' 'hr_swtpca' 'hr_cwtpca' 'hr_cwtpica' 'hr_cwtsobi' 'hr_ssa' 'hr_emd' ...
    'hr_eemd' 'hr_ceemd' 'hr_nmfcKL' 'hr_nmfcL2' 'hr_nmfcSparse'} ...
    ; ['RMSE' num2cell(heart_rmse2)] ...
    ; ['MAE' num2cell(heart_mae2)] ...
    ; {'','','','','','','','','','','','','',''}  ...
    ; {'Heart_Schmidt' 'hr_raw' 'hr_freqfilt' 'hr_swtpca' 'hr_cwtpca' 'hr_cwtpica' 'hr_cwtsobi' 'hr_ssa' 'hr_emd' ...
    'hr_eemd' 'hr_ceemd' 'hr_nmfcKL' 'hr_nmfcL2' 'hr_nmfcSparse'} ...
    ; ['RMSE' num2cell(heart_rmse)] ...
    ; ['MAE' num2cell(heart_mae)] ...
    ; {'','','','','','','','','','','','','',''} ...
    ; {'Lung' 'br_raw' 'br_freqfilt' 'br_swtpca' 'br_cwtpca' 'br_cwtpica' 'br_cwtsobi' 'br_ssa' 'br_emd' ...
    'br_eemd' 'br_ceemd' 'br_nmfcKL' 'br_nmfcL2' 'br_nmfcSparse'} ...
    ; ['RMSE' num2cell(lung_rmse)] ...
    ; ['MAE' num2cell(lung_mae)] ...
    ; {'','','','','','','','','','','','','',''} ] ...
    ,'VitalSignsAnalysis.xlsx','Sheet', 'Summary_per_min');

writecell([{'Heart_Springer' 'hr_raw' 'hr_freqfilt' 'hr_swtpca' 'hr_cwtpca' 'hr_cwtpica' 'hr_cwtsobi' 'hr_ssa' 'hr_emd' ...
    'hr_eemd' 'hr_ceemd' 'hr_nmfcKL' 'hr_nmfcL2'  'hr_nmfcSparse'} ...
    ; ['RMSE' num2cell(heart_rmse2/6)] ...
    ; ['MAE' num2cell(heart_mae2/6)] ...
    ; {'','','','','','','','','','','','','',''}  ...
    ; {'Heart_Schmidt' 'hr_raw' 'hr_freqfilt' 'hr_swtpca' 'hr_cwtpca' 'hr_cwtpica' 'hr_cwtsobi' 'hr_ssa' 'hr_emd' ...
    'hr_eemd' 'hr_ceemd' 'hr_nmfcKL' 'hr_nmfcL2' 'hr_nmfcSparse'} ...
    ; ['RMSE' num2cell(heart_rmse/6)] ...
    ; ['MAE' num2cell(heart_mae/6)] ...
    ; {'','','','','','','','','','','','','',''} ...
    ; {'Lung' 'br_raw' 'br_freqfilt' 'br_swtpca' 'br_cwtpca' 'br_cwtpica' 'br_cwtsobi' 'br_ssa' 'br_emd' ...
    'br_eemd' 'br_ceemd' 'br_nmfcKL' 'br_nmfcL2'   'br_nmfcSparse'} ...
    ; ['RMSE' num2cell(lung_rmse/6)] ...
    ; ['MAE' num2cell(lung_mae/6)] ...
    ; {'','','','','','','','','','','','','',''} ] ...
    ,'VitalSignsAnalysis.xlsx','Sheet', 'Summary_per_10sec');

fig=uifigure;
uitable(fig,'Data',readtable('VitalSignsAnalysis.xlsx','Sheet', 'Summary_per_10sec'))

