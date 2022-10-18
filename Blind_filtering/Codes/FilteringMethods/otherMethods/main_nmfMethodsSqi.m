% Initializing
clear
% close all
clc


%% Defining address of the folder containing the audio records
DA='..\..\..\chestRecords\currentDataset\';
dataSheet=readcell([DA '00dataSheet.xlsx']);
filesName=dataSheet([false; cell2mat(dataSheet(2:end,5))==1],1);

methodsName=["nmfc"]; 

TFs=["STFT","Q-transform"];
nmfTypes=["nmf_kl","nmf_general","nmf_sparse", "nmf_is_sparse", "nmf_ed_sparse5"];
reconstMethods=["Synthesis","No Mask","Filtering","Best Mask"];

%%
HeartSqiAllMethods=table;
LungSqiAllMethods=table;
TimeAllMethods=table;

for mn=methodsName
    method=struct;
    method.name=mn;
    if mn=="emd"
        for i=emdTypes
            for j=heartLocMethods
                method.emdType=i;
                method.heartLoc=j;
                [heartSqi, lungSqi, T] = getHeartLungSqiOfMethods(DA, filesName, method);
                TimeAllMethods = [TimeAllMethods table(T,'VariableNames', string([char(mn) '_' char(i) '_' char(j)]))];
                HeartSqiAllMethods=[HeartSqiAllMethods table(heartSqi,'VariableNames', string([char(mn) '_' char(i) '_' char(j)]))];
                LungSqiAllMethods = [LungSqiAllMethods table(lungSqi,'VariableNames', string([char(mn) '_' char(i) '_' char(j)]))];
            end
        end
    elseif mn=="nmfc"
       for i=TFs
           method.TF=i;
            for j=nmfTypes
                switch j
                    case "nmf_kl"
                        method.options_nmf.cf = 'kl';
                        method.options_nmf.sparsity = 0;
                    case "nmf_general"
                        method.options_nmf.cf = 'ed';
                        method.options_nmf.sparsity = 0;
                    case "nmf_sparse"
                        method.options_nmf.cf = 'kl';
                        method.options_nmf.sparsity = 2;
                   case "nmf_is_sparse"
                        method.options_nmf.cf = 'is';
                        method.options_nmf.sparsity = 2;
                    case "nmf_ed_sparse5"
                        method.options_nmf.cf = 'ed';
                        method.options_nmf.sparsity = 5;
                end
                for k=reconstMethods
                    method.reconst=k;
                    [heartSqi, lungSqi,T] = getHeartLungSqiOfMethods(DA, filesName, method);
                    TimeAllMethods=[TimeAllMethods table(T,'VariableNames', string([char(mn) '_' char(i) '_' char(j) '_' char(k)]))];
                    HeartSqiAllMethods=[HeartSqiAllMethods table(heartSqi,'VariableNames', string([char(mn) '_' char(i) '_' char(j) '_' char(k)]))];
                    LungSqiAllMethods = [LungSqiAllMethods table(lungSqi,'VariableNames', string([char(mn) '_' char(i) '_' char(j) '_' char(k)]))];
                end
            end
       end
    elseif mn == "ssa"
        [heartSqi, lungSqi,T] = getHeartLungSqiOfMethods(DA, filesName, method);
        heartSqi(heartSqi>5)=5; heartSqi(isnan(heartSqi))=1;
        lungSqi(lungSqi>5)=5; lungSqi(isnan(lungSqi))=1;
        TimeAllMethods = [TimeAllMethods table(T,'VariableNames', mn)];
        HeartSqiAllMethods = [HeartSqiAllMethods table(heartSqi,'VariableNames', mn)];
        LungSqiAllMethods = [LungSqiAllMethods table(lungSqi,'VariableNames', mn)];      
    elseif mn == "swtpca"
        method.N=6;
        method.wname='Bior1.3';
        Mdl=load ('MdlKnn.mat'); config = load ('config.mat');
        method.Mdl = Mdl.Mdl; method.config = config.config;
        [heartSqi, lungSqi,T] = getHeartLungSqiOfMethods(DA, filesName, method);
        TimeAllMethods = [TimeAllMethods table(T,'VariableNames', mn)];
        HeartSqiAllMethods = [HeartSqiAllMethods table(heartSqi,'VariableNames', mn)];
        LungSqiAllMethods = [LungSqiAllMethods table(lungSqi,'VariableNames', mn)];      
    elseif mn == "cwtpca"
        method.vpo = 4;
        method.wname = 'amor';
        method.fl = 50; method.fh=900;
        Mdl=load ('MdlKnn.mat'); config = load ('config.mat');
        method.Mdl = Mdl.Mdl; method.config = config.config;
        [heartSqi, lungSqi,T] = getHeartLungSqiOfMethods(DA, filesName, method);
        TimeAllMethods = [TimeAllMethods table(T,'VariableNames', mn)];
        HeartSqiAllMethods = [HeartSqiAllMethods table(heartSqi,'VariableNames', mn)];
        LungSqiAllMethods = [LungSqiAllMethods table(lungSqi,'VariableNames', mn)];      
    elseif mn == "cwtpica"
        method.vpo = 4;
        method.wname = 'amor'; 
        method.fl = 50; method.fh=900;
        method.AHR = 2.2;
        Mdl=load ('MdlKnn.mat'); config = load ('config.mat');
        method.Mdl = Mdl.Mdl; method.config = config.config;
        [heartSqi, lungSqi,T] = getHeartLungSqiOfMethods(DA, filesName, method);
        TimeAllMethods = [TimeAllMethods table(T,'VariableNames', mn)];
        HeartSqiAllMethods = [HeartSqiAllMethods table(heartSqi,'VariableNames', mn)];
        LungSqiAllMethods = [LungSqiAllMethods table(lungSqi,'VariableNames', mn)];      
    elseif mn == "cwtsobi"
        method.vpo = 4;
        method.wname = 'amor';
        method.fl = 50; method.fh=900;       
        method.tau = 'auto';
        Mdl=load ('MdlKnn.mat'); config = load ('config.mat');
        method.Mdl = Mdl.Mdl; method.config = config.config;
        [heartSqi, lungSqi,T] = getHeartLungSqiOfMethods(DA, filesName, method);
        TimeAllMethods = [TimeAllMethods table(T,'VariableNames', mn)];
        HeartSqiAllMethods = [HeartSqiAllMethods table(heartSqi,'VariableNames', mn)];
        LungSqiAllMethods = [LungSqiAllMethods table(lungSqi,'VariableNames', mn)]; 
    elseif mn == "freqfilt"   
        fs=4000;
        Hfr=[50 250];
        Lfr=[100 1000];
        [method.B_H,method.A_H] = butter(5,[2*Hfr(1)/fs 2*Hfr(2)/fs],'bandpass');
        [method.B_L,method.A_L] = butter(5,[2*Lfr(1)/fs 2*Lfr(2)/fs],'bandpass');
        [heartSqi, lungSqi,T] = getHeartLungSqiOfMethods(DA, filesName, method);
        TimeAllMethods = [TimeAllMethods table(T,'VariableNames', mn)];
        HeartSqiAllMethods = [HeartSqiAllMethods table(heartSqi,'VariableNames', mn)];
        LungSqiAllMethods = [LungSqiAllMethods table(lungSqi,'VariableNames', mn)]; 
    elseif mn == "raw"
        [heartSqi, lungSqi,T] = getHeartLungSqiOfMethods(DA, filesName, method);
        TimeAllMethods = [TimeAllMethods table(T,'VariableNames', mn)];
        HeartSqiAllMethods = [HeartSqiAllMethods table(heartSqi,'VariableNames', mn)];
        LungSqiAllMethods = [LungSqiAllMethods table(lungSqi,'VariableNames', mn)];         
    end
save temp1.mat   
end
writetable(TimeAllMethods, 'TimeAllMethods.xlsx');
writetable(HeartSqiAllMethods, 'HeartSqiNmfMethods.xlsx');
writetable(LungSqiAllMethods, 'LungSqiNmfMethods.xlsx');
% 
% M=array2table(mean(table2array(HeartSqiAllMethods)), 'VariableNames', HeartSqiAllMethods.Properties.VariableNames);
% Md=array2table(median(table2array(HeartSqiAllMethods)), 'VariableNames', HeartSqiAllMethods.Properties.VariableNames);


%%
%%%%% functions 

function [heartSqi, lungSqi, T] = getHeartLungSqiOfMethods(DA, filesName, method)
T=zeros(size(filesName,1),1);
for i= [] %1:size(filesName,1)

    %%% Display the progress percentage 
    disp([num2str(i) '/' num2str(size(filesName,1)) ', ' num2str(floor(100*i/size(filesName,1))) ' %']) 

    %%% Loading the signal
    [data,fs]=audioread([DA  '\' filesName{i}]);    %% read the record
    
    %%% Preprocessing - Downsampling to 4000 Hz
    data=downsample(data,fs./4000); fs=4000;
   
    %%% Preprocessing - normalizing
    data=normalize(data);
    
    %%% setting length to 10 sec
    data=data(1:10*fs);
    tic
    switch method.name
        case 'emd'
            [heartSound,lungSound]=emd_separation(data,fs, method.emdType, method.heartLoc);
            T(i) = toc;
        case 'nmfc'
            xhat= nmf_cluster_modified(data,fs,[], method.reconst, method.TF, [], method.options_nmf);         
            T(i) = toc;
            heartSound=xhat(:,1); lungSound=xhat(:,2);
        case 'ssa'
            [heartSound,lungSound]=singular_spectrum_analysis(data,fs);
            T(i) = toc;
        case 'cwtpca'
            S =cwtpca(data, fs, method.vpo, method.wname, method.fl, method.fh,  method.Mdl,  method.config);  
            T(i) = toc;
            heartSound=S{1}; lungSound=S{2};
        case 'swtpca'
            S = swtpca(data, fs,  method.N,  method.wname,(2:method.N),  method.Mdl,  method.config);
            T(i) = toc;
            heartSound=S{1}; lungSound=S{2};
        case 'cwtpica'
            S = cwtpica(data, fs,  method.vpo,  method.wname,  method.fl,  method.fh,  method.AHR,  method.Mdl,  method.config);          
            T(i) = toc;
            heartSound=S{1}; lungSound=S{2};
        case 'cwtsobi'
            S = cwtsobi(data, fs,  method.vpo,  method.wname,  method.fl,  method.fh,  method.tau,  method.Mdl,  method.config);      
            T(i) = toc;
            heartSound=S{1}; lungSound=S{2};
        case 'raw'
            heartSound = data(:);
            lungSound = data(:);
            T(i) = 0;
        case 'freqfilt'
            heartSound = filtfilt(method.B_H,method.A_H,data);
            lungSound = filtfilt(method.B_L,method.A_L,data);
            T(i) = toc;
    end
    
    
    
    
    
    heartSound=heartSound(:)./(6.*(std(heartSound(:),0,1))); % normalization
    lungSound=lungSound(:)./(15.*(std(lungSound(:),0,1))); % normalization
   
    try
        heartfeat(i,:)=get_all_SQIs_modified(heartSound, fs);
    catch
        warning('on');
        warning(['Heart feature extraction of case no.' num2str(i) ' is failed!'])
    end
    try
        lungfeat(i,:)=get_all_SQIs_modified(lungSound, fs);
    catch
        warning('on');
        warning(['Lung feature extraction of case no.' num2str(i) ' is failed!'])
    end
    save temp0.mat
end

[~,I]=sort(heartfeat.Properties.VariableNames);
ftrs=heartfeat(:,I); 
ftrs=table2array(ftrs);
load('heartAllFeatCVMdl.mat');
Sqi = predict(heartCVMdl,normalize(ftrs(:,heartFeatSelectIdx)')');
[~,I]=sort(filesName); 
heartSqi=Sqi(I);

[~,I]=sort(lungfeat.Properties.VariableNames);
ftrs=lungfeat(:,I); 
ftrs=table2array(ftrs);
load('lungAllFeatCVMdl.mat');
Sqi = predict(lungCVMdl,normalize(ftrs(:,lungFeatSelectIdx)')');
[~,I]=sort(filesName); 
lungSqi=Sqi(I);

 end





