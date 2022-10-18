%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%% running order of the main codes
%%% on a typical personal computers, the run may take about 48 hours.
clear 
close all
clc


main_TrainClassifiers
main_CompareClassifiers

main_FreqFiltButter
MainSWTPCA_autoCompClassif
MainCWTPCA_autoCompCalssif
MainCWTPiCA_autoCompClassif
MainCWTSOBI_autoCompClassif
main_nmfMethodsSqi
mainOtherMethods_filterAndSave

main_featureExtraction
main_featRegMdlTrain
main_getSqi

main_vitalSignError
main_sqiViolinGraph
main_samplesignals