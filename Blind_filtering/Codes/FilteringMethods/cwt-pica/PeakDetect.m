function Peaks=PeakDetect(x,ahr)



N=length(x);
Crng=floor(1./ahr);
Prng=floor(1.4*Crng);
Mrng=floor(.7*Crng);

[~,loc]=max(x(1:Prng));
locs=loc;
while loc+Prng<=N
    [~,I]=max(x(loc+Mrng:loc+Prng));
    loc=loc+Mrng+I-1;
    locs=[locs loc];
end
Peaks=false(N,1);
Peaks(locs)=true;


