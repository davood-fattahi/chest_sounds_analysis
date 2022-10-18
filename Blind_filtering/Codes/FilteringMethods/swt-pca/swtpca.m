function [S, varargout] =swtpca(data, fs, N, wname, decin, varargin)

% Single Channel Blind Source Separation Using SWTPCA method.
% 
% minimal syntax:
% [S,__]=swtpca(data, fs, N, wname, decin, p)
% [S,__]=swtpca(data, fs, N, wname, decin, Mdl, config)

% 
% Compelet syntax with all optional outputs:
% [S, Fltrs, frange, Us, U, p, Cmpnts, WCoefs, FrFltBnk] =swtpca(data, fs, N, wname, decin, p)
% [S, Fltrs, frange, Us, U, p, Cmpnts, WCoefs, FrFltBnk] =swtpca(data, fs, N, wname, decin, Mdl, config)
% 
% Inputs:
% data: tha input single channel data.
% fs: sampling frequency in Hz.
% N: no. of levels in SWT.
% wname: name of the utilized wavelet
% decin: index of decompositions which are involved in. e.g. [2 3 4].
% p: the desired components after applying the CA method. 
%   It can be a vector indexing the desired components for extracting one source. 
%   Or, can be in cell format whose each cell contains the components' index of each source.
% Mdl: a discriminant model provided for 'predict' function.
% config: the configuration of component classification. 
% 
% Outputs:
% p: index of independent sources, integers between 1 and size(X,1).
%   For one extracted source in output, it could be in a vector form.
%   For more than one extraxted sources, p is in cell format whose each cell 
%   contains correspounding index vector.
% 	Note: each 'extracted source' could be estimated from one or few
% 	'independent sources' indexed in p.
% S: the extracted sources, if p is vector (means one source is considered)
%   it will be a signal in vector form. But if p is in cell format (means more
%   than one sources are considered), S will be in cell format, too, whose 
%   each cell is correpounding to one source.
% Us: the vector of coefficients correspounding to each source; if p is a 
%   vector (means one source is considered) it will be in vector form. But if 
%   p is in cell format (means more than one sources are considered), Us will
%   also be in cell format, whose each cell is correpounding to one source.
% Filtrs: the estimated frequency filters correspounding to each source; if p is a 
%   vector (means one source is considered) it will be in vector form. But if 
%   p is in cell format (means more than one sources are considered), Fltrs will
%   be also in cell format, whose each cell is correpounding to one source.
% frange: frequancy range
% Cmpnts: the components after applying the CA method.
% Coefs: the wavelet coefficients.
% FrFltBnk: wavelet filterbank in frequancy domain.
% 
% 
% developed by Davood Fattahi (fattahi.d@gmail.com), 14/12/2020. 
% As a part of chest sound analysis project, under supervision of Faezeh MarzbanRad and Reza Sameni. 


%% Preprocessing - normalizing
data=normalize(data);

%%% freq. filters
    [FiltApprox, FiltDetail]=wavelet_filters(wname,N,fs); % discrete wavelet frequency filters
    FrFltBnk=[FiltDetail; FiltApprox(end,:)];
    FrFltBnk=FrFltBnk(decin,:);

%%% SWT
data_e1=data(1:end-mod(size(data,1),2^N));  % correcting the length of signal
SWC=swt(data_e1,N,wname); % SWT
%%    
%%% PCA
WCoefs=SWC(decin,:);
[~,U,Cmpnts,~]=mypca(WCoefs,size(WCoefs,1));

%% Providing other outputs
if nargin==6
    p=varargin{1};
elseif nargin==7
    [lbls, nh, nl]=componentClassifyTest(Cmpnts', fs, varargin{1}, varargin{2});
    p=cell(2,1);
    p{1}=find(contains(lbls,'h'));
    p{2}=find(contains(lbls,'l'));
end

% if p is empty
if isempty(p)
    p=1;
end
lp=length(p);

if iscell(p)
    S=cell(lp,1);
    Us=cell(lp,1);
    Fltrs=cell(lp,1);
    for i=1:lp
        SWCr=zeros(size(SWC)); 
        SWCr(decin,:)=U(:,p{i})*Cmpnts(p{i},:);
        S{i}=iswt(SWCr,wname);
        Us{i}=U(:,p{i});
        
        dd=zeros(length(data_e1),1); dd(ceil(length(data_e1)/2))=1; % impulse
        swc=swt(dd,N,wname); % SWT
        swcr=zeros(size(swc));
        swcr(decin,:)=U(:,p{i})*U(:,p{i})'*swc(decin,:);
        filters=abs(fft(iswt( swcr ,wname),ceil(0.5*length(data)))); % impulse response
        Fltrs{i}=filters(1:ceil(0.5*length(filters))); % one-side frequencies

%         Fltrs{i}=abs(sum((U(:,p{i})*U(:,p{i})'*FrFltBnk)));
    end
else
    SWCr=zeros(size(SWC)); 
    SWCr(decin,:)=U(:,p)*Cmpnts(p,:);
    S=iswt(SWCr,wname);
    Us=U(:,p);

    dd=zeros(length(data_e1),1); dd(ceil(length(data_e1)/2))=1; % impulse
    swc=swt(dd,N,wname); % SWT
    swcr=zeros(size(swc));
    swcr(decin,:)=U(:,p)*U(:,p)'*swc(decin,:);
    filters=abs(fft(iswt( swcr ,wname),ceil(0.5*length(data)))); % impulse response
    Fltrs=filters(1:ceil(0.5*length(filters))); % one-side frequencies

%     Fltrs=abs(sum((U(:,p)*U(:,p)'*FrFltBnk)));
end
frange = .5*fs*(1:(ceil(0.5*length(filters))))/(ceil(0.5*length(filters)));



if nargout>=2
varargout{1}=Fltrs;
end
if nargout>=3
varargout{2}=frange;
end
if nargout>=4
varargout{3}=Us;
end
if nargout>=5
varargout{4}=U;
end
if nargout>=6
varargout{5}=p;
end
if nargout>=7
varargout{6}=Cmpnts;
end
if nargout >=8
varargout{7}=WCoefs;
end
if nargout>=9
varargout{8}=FrFltBnk;
end
if nargout>=10
varargout{9}=nh;
end
if nargout>=11
varargout{10}=nl;
end

end






