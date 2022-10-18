function [U,S,Xh]=pwsobi(X,nw,nvrlp,tau,p,varargin)

% piecewise second order blind identification
% 
% Inputs:
% X: the input signal, (T-F decomposed or multi-channel)
% nw: moving window size.
%   --> Scalar for equi-size windows; 
%   --> Vector for variable size windows. For example: nw=[200 300] means two windows with 200 and
%       300 samples size (overlap is not allowed in this form);
%   --> Two columns matrix, 1st column for start point and 2-nd column for
%       end point of each segment.
% nvrlp:overlapping size
% tau: time-lags in SOBI. 
%   --> scalar: number of one sample time-lags (or number of cov. matrices 
%       to be diagonalized in SOBI. For example, 2000 means 2000 time-lags 
%       from 1 to 2000 samples.
%   --> vector: vector of time-lags (samples) in SOBI. For example, 
%       (2;200;2000) means 3 time lags of 2, 200 and 2000 samples.
% p: index of independent sources, integers between 1 and size(X,1).
%   For one extracted source in output, it could be in a vector form.
%   For more than one extraxted sources, p is in cell format and each cell 
%   contains correspounding index vector.
% 	Note: each 'extracted source' could be estimated from one or few
% 	'independent sources' indexed in p.
% vp: property of signal values used for cov matrix calculation; 
%       'defualt': signal values
%       'abs': absolute values af signal
%       'real': real values of signal
%       'imag': imaginary vaues of signal
%       'renv': envelope of the real signal
% 
% 
% Outputs:
% U: the 3-D matrix of decomposing coeficients, the 1th dimension (each row)
%   for each channel (or frequency band in T-F decomposition), the 2nd 
%   dimension (each column) for each independent source, and the 3rd dimension for each 
%   time sample.
% S: the extracted sources.
% Xh: the reconstructed signal using the sources.
% Note: If the extracted sources are more than two, U will be in
%   cell format.
% 
% 
%
% developed by Davood Fattahi (fattahi.d@gmail.com), 18/07/2020, based on 
% the SOBI algorithm developed by A. Belouchrani and A. Cichocki
% As a part of chest sound analysis project, under supervision of Faezeh MarzbanRad and Reza Sameni. 
%%
d=size(X,1); n=size(X,2); lp=length(p); % just for convenience
% X=X-repmat(mean(X,2),1,n); % make the signal zero-mean 


%% property of signal values used for cov matrix calculation
if isempty(varargin)||strcmp(varargin{1},'default')
    Y=X;
elseif strcmp(varargin{1},'abs')
    Y=abs(X);
elseif strcmp(varargin{1},'real')
    Y=real(X);
elseif strcmp(varargin{1},'imag')
    Y=imag(X);
elseif strcmp(varargin{1},'renv')
    Y=abs(hilbert(real(X)));
else
    Y=X;    
end

%% if p is empty
if isempty(p)
    p=1;
end

%% if tau is empty or out of range
if isempty(tau)
    tau=100;
end
  









%%% windowing
if ~isvector(nw)
    wi=nw;
else
    nw=nw(:); % size of the sliding windows
    if length(nw)==1
        K=ceil((n-nw+1)./(nw-nvrlp));  % total no. of windows

        wi=zeros(K,2); % the windows start and end index 
        for k=1:K
            if (k-1)*(nw-nvrlp)+nw <= n
            wi(k,:)=[(k-1)*(nw-nvrlp)+1 (k-1)*(nw-nvrlp)+nw];
            else
            wi(k,:)=[(k-1)*(nw-nvrlp)+1 n];    
            end    
        end

    elseif isempty(nw) % if nw is empty, nw will be the size of each beats
        wi=find(Peaks); wi=wi(:); wi=[1;wi;length(Peaks)]; wi=[wi(1:end-1) wi(2:end)-1];

    elseif length(nw)>1 % if nw is vector of windows sizes
        nw=[1; nw];
        for k=1:length(nw)-1
            wi(k,:)=[sum(nw(1:k)) sum(nw(1:k+1))-1];
        end
    end 
end
wi=wi(wi(:,2)+min(tau)<n,:);
%%% pre allocation
if iscell(p)
    for i=1:lp
        S{i}=zeros(length(p{i}),n);
        Xh{i}=zeros(d,n);
        U{i}=zeros(d,length(p{i}),n);
    end
else
    S=zeros(lp,n);
    Xh=zeros(d,n);
    U=zeros(d,lp,n);
end
rp=zeros(1,n);
K=size(wi,1);
Prcnt=0;

for k=1:K % for each sliding window
    %%% display the percentage of progressing
    prcnt=floor(100*k/K);
    if prcnt > Prcnt
        Prcnt=prcnt;
        disp([num2str(k) '/' num2str(K) ', ' num2str(prcnt) ' %']) 
    end
    y=Y(:,wi(k,1):wi(k,2)); % signal window
    
    %%% SOBI
    [u_d,s,~]=wsobi(Y,d,tau(tau+wi(k,2)<=n), wi(k,:));
    
        %%% selecting p eig vectors
    if ~iscell(p)  
        u_p=u_d(:,p(:));
        U(:,:,wi(k,1):wi(k,2))=U(:,:,wi(k,1):wi(k,2))+repmat(u_p,1,1,wi(k,2)-wi(k,1)+1);
        S(:,wi(k,1):wi(k,2))=S(:,wi(k,1):wi(k,2))+s(p(:),:);
        xh=u_p*s(p(:),:);
        Xh(:,wi(k,1):wi(k,2))=Xh(:,wi(k,1):wi(k,2))+xh;
    else
        for i=1:lp
            u_p=u_d(:,p{i});
            U{i}(:,:,wi(k,1):wi(k,2))=U{i}(:,:,wi(k,1):wi(k,2))+repmat(u_p,1,1,wi(k,2)-wi(k,1)+1);
            S{i}(:,wi(k,1):wi(k,2))=S{i}(:,wi(k,1):wi(k,2))+s(p{i},:);
            xh=u_p*s(p{i},:);
            Xh{i}(:,wi(k,1):wi(k,2))=Xh{i}(:,wi(k,1):wi(k,2))+xh;
        end    
     end
    
    rp(wi(k,1):wi(k,2))=rp(wi(k,1):wi(k,2))+1; % times of repeating
end
rp(rp==0)=1; % assure no zero values (does not happen usually)

%%% averaging over overlaped samples
if iscell(p)
    for i=1:lp
        S{i}=S{i}./repmat(rp,length(p{i}),1);
        U{i}=U{i}./(repmat(permute(rp,[1 3 2]),d,length(p{i}),1));
        Xh{i}=Xh{i}./repmat(rp,d,1);
    end
else
    S=S./repmat(rp,lp,1);
    U=U./(repmat(permute(rp,[1 3 2]),d,lp,1));
    Xh=Xh./repmat(rp,d,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    


    
    
