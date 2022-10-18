function [allFeatures]=collectAllFeatures(DA,filesName)


if ischar(DA)
    dacell=cell(size(filesName));
    dacell(:)={DA};
    DA=dacell;
end
allFeatures=table;
for i=1:size(filesName,1)
    lineLength=fprintf('file name: %s, file number: %d \n', filesName{i} ,i);
    % load the audio file
    if ~(DA{i}(end)=='\' || DA{i}(end)=='/')
        DA{i}(end+1)='\';
    end
    [audio_data, Fs]=audioread([DA{i} filesName{i}]);

    % extract the corresonding features
    warning('on');
    try
        features=get_all_SQIs_modified(audio_data, Fs);
        allFeatures(i,:) = [cell2table(filesName(i)) features];
    catch
        warning(['Error in ' filesName{i} ', the file is skipped!']);
        allFeatures(i,1) = cell2table(filesName(i));
    end
end


