

function [classifyResult]= challenge(recordName)

load('Springer_B_matrix.mat');
load('Springer_pi_vector.mat');
load('Springer_total_obs_distribution.mat');
[data, Fs] = audioread([recordName, '.wav']);

audio_data2 =butterworth_high_pass_filter(data,2,700,Fs);
sigma=std(audio_data2);
if 4*sigma >= 0.17
    classifyResult=0;
else
    %% Load data and resample data
    springer_options   = default_Springer_HSMM_options;
    % load data
    [PCG, Fs1] = audioread([recordName  '.wav']);
    % resample to springer_options.audio_Fs (1000 Hz)
    PCG_resampled  = resample(PCG,springer_options.audio_Fs,Fs1);
    
    %% Running runSpringerSegmentationAlgorithm.m to obtain the assigned_states
    [assigned_states,denoised_signal] = runSpringerSegmentationAlgorithm2(PCG_resampled , springer_options.audio_Fs, Springer_B_matrix, Springer_pi_vector, Springer_total_obs_distribution, false); % obtain the locations for S1, systole, s2 and diastole
    [signal_s1,signal_s2, signal_systole,signal_diastole]=seperate_states_2(assigned_states,denoised_signal,false);
    
    [features,A,a]  = extractFeaturesFromHsIntervals2(assigned_states,PCG_resampled,signal_systole,signal_diastole,signal_s1,signal_s2);
    
    %% extracting  11 wt entropy features
    %article: christer ahlstrom
    [entropy_matrix]=mean_beat_wt_entropy(denoised_signal,a);
    
    %% extracting 9 wt features
    %article:christer ahlstrom
    [sum_matrix_mean]=absolute_sum(denoised_signal,a);
    
    %% extracting 9 shannon featueres
    [shannon_matrix]=shannon_en(denoised_signal,a);
    
    %% calculating mfcc features
    %  [MFCCs_s1  MFCCs_s2 MFCCs_systole MFCCs_diastole ]=calculating_mfcc(signal_s1 ,signal_s2 ,signal_systole, signal_diastole,A)
    [mean_mfcc,mean_mfcc_systole ,mean_mfcc_diastole]=calculating_mfcc(signal_s1 ,signal_s2 ,signal_systole, signal_diastole,a);
    feature_taki_mfcc=max( mean_mfcc_systole(4) ,mean_mfcc_diastole(4) );
    
    %[mfcc_s1,mfcc_s2,mfcc_sys,mfcc_dis]=seperate_mfcc(mean_mfcc);
    
    %% Running extractFeaturesFromHsIntervals.m to obtain the features for normal/abnormal heart sound classificaiton
    
    %[mean_entropy_systole,std_entropy_systole,skewness_entropy_systole,kurtosis_entropy_systole,mean_entropy_diastole,std_entropy_diastole,kurtosis_entropy_diastole,skewness_entropy_diastole,mean_power_sys,std_power_sys,kurtosis_power_sys,skewness_power_sys,mean_power_dias,std_power_dias,kurtosis_power_dias,skewness_power_dias]=power_entropy(a,signal_systole,signal_diastole)
    power_entropy_mean=power_entropy2(a,signal_systole,signal_diastole);
    
    all_features=[entropy_matrix,feature_taki_mfcc,mean_mfcc(:,1)',mean_mfcc(:,2)',mean_mfcc(:,3)',mean_mfcc(:,4)',power_entropy_mean',shannon_matrix,sum_matrix_mean];
end
end


function [features,A,a] = extractFeaturesFromHsIntervals2(assigned_states, PCG,signal_systole,signal_diastole,signal_s1,signal_s2)
% This function calculate 20 features based on the assigned_states by running "runSpringerSegmentationalgorithm.m" function
%
% INPUTS:
% assigned_states: the array of state values assigned to the sound recording.
% PCG: resampled sound recording with 1000 Hz.
%
% OUTPUTS:
% features: the obtained 20 features for the current sound recording
%
% Written by: Chengyu Liu, January 22 2016
%             chengyu.liu@emory.edu
%
% Last modified by:
%
% $$$$$$ IMPORTaNT
% Please note: the calculated 20 features are only some pilot features, some features maybe
% helpful for classifying normal/abnormal heart sounds, some maybe
% not. You need re-construct the features for a more accurate classification.


%% We just assume that the assigned_states cover at least 2 whole heart beat cycle
indx = find(abs(diff(assigned_states))>0); % find the locations with changed states

% for some recordings, there are state zeros at the beginning of assigned_states
if assigned_states(1)>0
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

indx2 = indx(K:end);
rem = mod(length(indx2),4);
indx2(end-rem+1:end) = [];
% a is N*4 matrix, the 4 columns save the beginnings of S1, systole, S2 and diastole in the same heart cycle respectively
A  = reshape(indx2,4,length(indx2)/4)';
[a]=mofify_A(A,signal_systole,signal_diastole,signal_s1,signal_s2);
if a(end,5)==0
    a(end,:)=[];
end
b=isempty(a);
if b==1
    features=0;
else
    %% Feature calculation
    % mean value of RR intervals
    m_RR        = round(mean((a(:,5)-a(:,1))));
    % standard deviation (SD) value of RR intervals
    sd_RR       = round(std(a(:,5)-a(:,1)));
    % mean value of S1 intervals
    mean_IntS1  = round(mean(a(:,2)-a(:,1)));
    % SD value of S1 intervals
    sd_IntS1    = round(std(a(:,2)-a(:,1)));
    % mean value of S2 intervals
    mean_IntS2  = round(mean(a(:,4)-a(:,3)));
    % SD value of S2 intervals
    sd_IntS2    = round(std(a(:,4)-a(:,3)));
    % mean value of systole intervals
    mean_IntSys = round(mean(a(:,3)-a(:,2)));
    % SD value of systole intervals
    sd_IntSys   = round(std(a(:,3)-a(:,2)));
    % mean value of diastole intervals
    mean_IntDia = round(mean(a(:,5)-a(:,4)));
    % SD value of diastole intervals
    sd_IntDia   = round(std((a(:,5)-a(:,4))));
    
    R_SysRR=zeros(size(a,1));
    R_DiaRR=zeros(size(a,1));
    R_SysDia=zeros(size(a,1));
    P_S1=zeros(size(a,1));
    P_Sys=zeros(size(a,1));
    P_S2=zeros(size(a,1));
    P_Dia=zeros(size(a,1));
    P_SysS1=zeros(size(a,1));
    P_DiaS2=zeros(size(a,1));
    for i=1:size(a,1)
        R_SysRR(i)  = (a(i,3)-a(i,2))/(a(i,5)-a(i,1))*100;
        R_DiaRR(i)  = (a(i,5)-a(i,4))/(a(i,5)-a(i,1))*100;
        R_SysDia(i) = R_SysRR(i)/R_DiaRR(i)*100;
        
        P_S1(i)     = sum(abs(PCG(a(i,1):a(i,2))))/(a(i,2)-a(i,1));
        P_Sys(i)    = sum(abs(PCG(a(i,2):a(i,3))))/(a(i,3)-a(i,2));
        P_S2(i)     = sum(abs(PCG(a(i,3):a(i,4))))/(a(i,4)-a(i,3));
        P_Dia(i)    = sum(abs(PCG(a(i,4):a(i,5))))/(a(i,5)-a(i,4));
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
    
    % mean value of the interval ratios between systole and RR in each heart beat
    m_Ratio_SysRR   = mean(R_SysRR);
    % SD value of the interval ratios between systole and RR in each heart beat
    sd_Ratio_SysRR  = std(R_SysRR);
    % mean value of the interval ratios between diastole and RR in each heart beat
    m_Ratio_DiaRR   = mean(R_DiaRR);
    % SD value of the interval ratios between diastole and RR in each heart beat
    sd_Ratio_DiaRR  = std(R_DiaRR);
    % mean value of the interval ratios between systole and diastole in each heart beat
    m_Ratio_SysDia  = mean(R_SysDia);
    % SD value of the interval ratios between systole and diastole in each heart beat
    sd_Ratio_SysDia = std(R_SysDia);
    
    % avoid the flat line signal
    indx_sys = find(P_SysS1>0 & P_SysS1<100);
    if length(indx_sys)>1
        % mean value of the mean absolute amplitude ratios between systole period and S1 period in each heart beat
        m_amp_SysS1  = mean(P_SysS1(indx_sys));
        % SD value of the mean absolute amplitude ratios between systole period and S1 period in each heart beat
        sd_amp_SysS1 = std(P_SysS1(indx_sys));
    else
        m_amp_SysS1  = 0;
        sd_amp_SysS1 = 0;
    end
    indx_dia = find(P_DiaS2>0 & P_DiaS2<100);
    if length(indx_dia)>1
        % mean value of the mean absolute amplitude ratios between diastole period and S2 period in each heart beat
        m_amp_DiaS2  = mean(P_DiaS2(indx_dia));
        % SD value of the mean absolute amplitude ratios between diastole period and S2 period in each heart beat
        sd_amp_DiaS2 = std(P_DiaS2(indx_dia));
    else
        m_amp_DiaS2  = 0;
        sd_amp_DiaS2 = 0;
    end
    
    features = [m_RR sd_RR  mean_IntS1 sd_IntS1  mean_IntS2 sd_IntS2  mean_IntSys sd_IntSys  mean_IntDia sd_IntDia m_Ratio_SysRR sd_Ratio_SysRR m_Ratio_DiaRR sd_Ratio_DiaRR m_Ratio_SysDia sd_Ratio_SysDia m_amp_SysS1 sd_amp_SysS1 m_amp_DiaS2 sd_amp_DiaS2];
end
end

function [a]=mofify_A(A,signal_systole,signal_diastole,signal_s1,signal_s2)
for i=1:size(A,1)-1
    A(i,5)=A(i+1,1);
end
%%
indicies_p_s_f=0;
indicies_p_d_f=0;

%% systole
j=1;
seperate_systoles=cell(size(A,1),1);
for i=1:size(A,1)
    seperate_systoles{j,1}=signal_systole(A(i,2):A(i,3));
    j=j+1;
end

m=1;
power_sys=0;
for k=1:size(A,1)
    power_sys(m)=bandpower(seperate_systoles{k,1});
    m=m+1;
end

% remove outliers for power_sys
j=1;
z=zscore(power_sys);
indicies_p_s=find(abs(z)>1);
for i=1:length(indicies_p_s)
    if power_sys(indicies_p_s(i))>=mean(power_sys)
        indicies_p_s_f(j)=indicies_p_s(i);
        j=j+1;
    end
end
%% diastole
j=1;
seperate_diastoles=cell(size(A,1),1);
for i=1:size(A,1)
    seperate_diastoles{j,1}=signal_diastole(A(i,4):A(i,5));
    j=j+1;
end
m=1;
power_dias=0;
for k=1:size(A,1)
    power_dias(m)=bandpower(seperate_diastoles{k,1});
    m=m+1;
end

% remove outliers for power_dias
j=1;
z=zscore(power_dias);
indicies_p_d=find(abs(z)>1);
for i=1:length(indicies_p_d)
    if power_dias(indicies_p_d(i))>=mean(power_dias)
        indicies_p_d_f(j)=indicies_p_d(i);
        j=j+1;
    end
end
%% finding bad s1
j=1;
seperate_s1=cell(size(A,1),1);

for i=1:size(A,1)
    seperate_s1{j,1}=signal_s1(A(i,1):A(i,2));
    j=j+1;
end

m=1;
power_s1=0;
for k=1:size(A,1)
    power_s1(m)=bandpower(seperate_s1{k,1});
    m=m+1;
end
indicies_p_s1_f=0;
j=1;

% remove outliers for power_s1
z=zscore(power_s1);
indicies_p_s1=find(abs(z)>1);
for i=1:length(indicies_p_s1)
    if power_s1(indicies_p_s1(i))<=mean(power_s1)
        indicies_p_s1_f(j)=indicies_p_s1(i);
        j=j+1;
    end
end
%% finding bad s2
j=1;
seperate_s2=cell(size(A,1),1);

for i=1:size(A,1)
    seperate_s2{j,1}=signal_s2(A(i,3):A(i,4));
    j=j+1;
    
end
m=1;
power_s2=0;
for k=1:size(A,1)
    power_s2(m)=bandpower(seperate_s2{k,1});
    m=m+1;
end

j=1;
indicies_p_s2_f=0;

% remove outliers for power_s2
z=zscore(power_s2);
indicies_p_s2=find(abs(z)>1);
for i=1:length(indicies_p_s2)
    if power_s2(indicies_p_s2(i))<=mean(power_s2)
        indicies_p_s2_f(j)=indicies_p_s2(i);
        j=j+1;
    end
end

%%

bad_cycles=[indicies_p_s_f, indicies_p_d_f,indicies_p_s1_f,indicies_p_s2_f];
bad_cycles=nonzeros(unique(sort(bad_cycles)))';

%% A_MODIFY
a=A;
a(bad_cycles,:)=[];
end





