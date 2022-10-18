Time=readtable('TimeAllMethods.xlsx');
meanTime = varfun(@mean, Time);
stdTime = varfun(@std, Time);

writetable(meanTime, 'meanStdTimeAllMethods.xlsx', 'Sheet', 'mean')
writetable(stdTime, 'meanStdTimeAllMethods.xlsx', 'Sheet', 'std')


