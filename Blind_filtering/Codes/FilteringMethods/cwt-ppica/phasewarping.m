function [y, ph]=phasewarping(x,Phase,bins)

% wraping the semiperiodic signal to the same phase
% 
% inputs:
% x: input semi-perodic signal
% phase: real value phase of the input signal, between [-pi pi]
% bins: the number of samples considered for each period of output
%
% outputs:
% y: the output phase-wraped signal.
%
% Daood Fatthi, fattahi.d@gmail.com, 10/03/2020
%%

P=find((Phase(1:end-1)-Phase(2:end))>0); P=[0 ;P(:); length(Phase)];
N=length(P);

ph=((2*pi)/bins:(2*pi)/bins:2*pi)-pi;

y=nan(bins,N-1);
for i=1:N-1
for j=1:bins % for each sample
    [~,a]=min(abs(ph(j)-Phase(P(i)+1:P(i+1))));
    y(j,i)=x(a+P(i));
end
end
