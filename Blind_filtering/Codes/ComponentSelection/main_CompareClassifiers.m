clear
% close all
clc

fs=4000;

load componentSelectionData.mat


[Mdl, config] = cmponentClassifyTrain(Cmpnnts,Lbls,fs);
% save('Mdl.mat','Mdl');
save('config.mat','config');


load config.mat
load 'MdlKnn.mat'
c=Cmpnnts(:,1:15);
labels=componentClassifyTest(c, fs, Mdl, config);

for i=1:size(Cmpnnts,2)
    tic;
    labels=componentClassifyTest(Cmpnnts(:,i), fs, Mdl, config);
    te(i)=toc;
end
M=mean(te(:)); S=std(te(:));
Tbl=table(M, S, 'variableName', {'Mean (sec)', 'STD (sec)'});
writetable(Tbl, 'KnnSelectorTimeCost.csv');

load config.mat
load 'MdlKnn.mat'
for i=1:10
    CVMdl = crossval(Mdl);
    classErrorKnn(i,:) = kfoldLoss(CVMdl,'Mode','individual' );
end

load 'MdlLda.mat'
for i=1:10
    CVMdl = crossval(Mdl);
    classErrorLda(i,:) = kfoldLoss(CVMdl,'Mode','individual' );
end

load 'MdlSvm.mat'
for i=1:10
    CVMdl = crossval(Mdl);
    classErrorSvm(i,:) = kfoldLoss(CVMdl,'Mode','individual' );
end

load 'MdlBys.mat'
for i=1:10
    CVMdl = crossval(Mdl);
    classErrorBys(i,:) = kfoldLoss(CVMdl,'Mode','individual' );
end

load 'MdlDecTree.mat'
for i=1:10
    CVMdl = crossval(Mdl);
    classErrorDecTree(i,:) = kfoldLoss(CVMdl,'Mode','individual' );
end
Method={'KNN'; 'LDA'; 'SVM'; 'Bys'; 'Tree'};
Mean=1-[mean(classErrorKnn(:)); mean(classErrorLda(:)); mean(classErrorSvm(:)); mean(classErrorBys(:)); mean(classErrorDecTree(:))];
Std=[std(classErrorKnn(:)); std(classErrorLda(:)); std(classErrorSvm(:)); std(classErrorBys(:)); std(classErrorDecTree(:))];

Tb2=table(Method, Mean, Std);
writetable(Tb2, 'CrossValAcc.csv');



fig1=uifigure;
uitable(fig1,'Data',TB1)


fig2=uifigure;
uitable(fig2,'Data',TB2)


