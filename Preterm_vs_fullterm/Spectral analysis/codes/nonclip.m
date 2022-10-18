function Y=nonclip(data,fs,wind,tr,Tr,piece)
% this function find the clipped segments in the input audio signal, and exclude them.


% Input:
% data: column signal
% fs: sampling frequency
% wind: time window
% tr: normalized threshold for detecting clipped samples 
% Tr: normalized threshold for detecting clipped time windows
% piece: type of selected pieces; can be 'longest' and 'concat'.
%     'longest': finds and select the longest nonclipped segment.
%     'concat': finds all the nonclipped segments and concatanates them
%        together.



%  Davood Fattahi 07/16/2019
%%%%%%%%%%%%%%%%%%%%%


wind =floor(wind*fs); % 1000 ms
Data=reshape(data(1:floor(size(data,1)/wind)*wind), wind,[]);

switch piece
    case 'longest'    
    k=1;
    s=[];
    ss=[];
    for j=1:size(Data,2)
        tR=tr*(max(Data(:,j))-min(Data(:,j)));
        mx=max(Data(:,j))-tR; mn=min(Data(:,j))+tR;
        [r, ~]=find(Data(:,j)>mx | Data(:,j)<mn);
        Clp=zeros(wind,1); Clp(r)=1;
        p(j) =sum(Clp,1)./wind;
        if p(j)<Tr
            s=[s ;j];
            if size(s,1)>k
                k=size(s,1);
                ss=s;
            end
        else
            s=[];
        end
    end
    
    case 'concat'
        for j=1:size(Data,2)
        tR=tr*(max(Data(:,j))-min(Data(:,j)));
        mx=max(Data(:,j))-tR; mn=min(Data(:,j))+tR;
        [r, ~]=find(Data(:,j)>mx | Data(:,j)<mn);
        Clp=zeros(wind,1); Clp(r)=1;
        p(j) =sum(Clp,1)./wind;
        end
        [~,ss]=find(p<Tr);
end
    Y=Data(:,ss);
    Y=reshape(Y,[],1);




    figure
    plot(data); 
    P=zeros(1,size(Data,2)); 
    P(ss)=max(Y); 
    P=repmat(P,wind,1);
    P=reshape(P,1,[]);
    hold on; plot(P,'LineWidth',2);
%     figure
%     plot(Y)






