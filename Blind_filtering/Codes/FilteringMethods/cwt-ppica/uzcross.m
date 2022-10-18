function a=uzcross(x)

% upward zerro crossing detection
% 
% Davood Fattahi, fattahi.d@gmail.com,  2019/12/10
%%
xm=x(1:end-1);
xp=x(2:end);
a=find((xp.*xm)<=0);
a=a((xm(a)-xp(a))<=0); % upward zerocrossing detection
a((a-circshift(a,1))==1)=[]; % removing repeated ones



