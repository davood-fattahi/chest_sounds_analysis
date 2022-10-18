function [T0,T1] = SynchBeatTimes(peaks)

% non warped beat synchronization ...

ip=find(peaks);
RR=(ip(2:end)-ip(1:end-1));
RR=[RR(1); RR; RR(end)];

T0=1:ip(end-1);
T1=T0;

for i=1:length(ip)
    ib= (T0>(ip(i)-ceil(.45*RR(i)))) & (T0<(ip(i)+ceil(.55*RR(i+1))));
    T1(ib)=T0(ib)+RR(i+1);
end