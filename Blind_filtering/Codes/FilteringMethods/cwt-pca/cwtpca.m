function [S, varargout] =cwtpca(data, fs, vpo, wname, fl, fh, varargin)
% Single Channel Blind Source Separation Using CWTPCA method.

% minimal syntax:
% [S,__]=cwtpca(data, fs, vpo, wname, fl, fh, p)
% [S,__]=cwtpca(data, fs, vpo, wname, fl, fh, Mdl, config)
% 
% 
% Compelet syntax with all optional outputs:
% [S, Fltrs, frange, Us, U, p, Cmpnts, WCoefs, FrFltBnk] =cwtpca(data, fs, vpo, wname, fl, fh, p)
% [S, Fltrs, frange, Us, U, p, Cmpnts, WCoefs, FrFltBnk] =cwtpca(data, fs, vpo, wname, fl, fh, Mdl, config)
% 
% Inputs:
% data: tha input single channel data.
% fs: sampling frequency in Hz.
% vpo: voices per octave in CWT.
% wname: name of the utilized wavelet.
% fl and fh: lower and upper frequency limits in CWT.
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


%% cwt
%%% filterbank generation
fb = cwtfilterbank('SignalLength',length(data),'Wavelet',wname, ...
'VoicesPerOctave',vpo, 'SamplingFrequency',fs, ...
'FrequencyLimits',[fl fh]); 

[FrFltBnk, ~]=freqz(fb); % cwt frequency filters

%%% cwt coefficients (complex value) 
[WCoefs, f]=cwt(data,'FilterBank',fb);


%%    
%%% PCA
[~,U,Cmpnts,~]=mypca(WCoefs,size(WCoefs,1));

%% Providing other outputs
if nargin==7
    p=varargin{1};
elseif nargin==8
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
        S{i}=icwt(U(:,p{i})*Cmpnts(p{i},:),wname,f,[f(end) f(1)]);
        Us{i}=U(:,p{i});
%         Fltrs{i}=abs(sum((U(:,p{i})*U(:,p{i})'* FrFltBnk))); %%  need to be normalized
%         Fltrs{i}=abs(fft(icwt((U(:,p{i})*U(:,p{i})'* wavelets(fb)),wname,f,[f(end) f(1)]),ceil(.5*length(data))));
        dd=zeros(length(data),1); dd(ceil(length(data)/2))=1; % impulse
        filters=abs(fft(icwt((U(:,p{i})*U(:,p{i})'* cwt(dd,'FilterBank',fb)),wname,f,[f(end) f(1)]),ceil(0.5*length(data)))); % impulse response
        Fltrs{i}=filters(1:ceil(0.5*length(filters))); % one-side frequencies
    end
elseif isvector(p)
    S=icwt(U(:,p)*Cmpnts(p,:),wname,f,[f(end) f(1)]);
    Us=U(:,p);
    dd=zeros(length(data),1); dd(ceil(length(data)/2))=1; % impulse
    filters=abs(fft(icwt((U(:,p{i})*U(:,p{i})'* cwt(dd,'FilterBank',fb)),wname,f,[f(end) f(1)]),ceil(0.5*length(data)))); % impulse response
    Fltrs=filters(1:ceil(0.5*length(filters))); % one-side frequencies
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






