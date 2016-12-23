clear;
EbNodB=6;
%EsNodB=10*log10(Rm*Rc£©+EbNodB
%EsNo=Rm*Rc*EbNo
scale=1000000;
x=randi([0,1],scale,1);
mod=2*x-1;

hChan = comm.AWGNChannel('EbNo',EbNodB);
y=step(hChan,mod);
d=y-mod;

y1=awgn(mod,EbNodB+3);
d1=y1-mod;

sigma2=1/2*10^(-EbNodB/10);
sigma=sqrt(sigma2);
noise=sigma*randn(scale,1);
y2=mod+noise;
d2=y2-mod;

figure(1);
hist(d,1000);
grid on;
figure(2);
hist(d1,1000);
grid on;
figure(3);
hist(d2,1000);
grid on;
% figure(3);
% a=-2:0.01:2;
% b=1/sqrt(2*pi)/sigma*exp(-a.^2/2/sigma^2);
% plot(a,b,'r');



xcomp=mod*exp(pi/4*1i);
y3=awgn(xcomp,EbNodB);
d3=abs(y3-xcomp);
% sigma2=sigma2*2;
% sigma=sqrt(sigma2);
noise4=(sigma*randn(scale,1)+1i*sigma*randn(scale,1));
y4=xcomp+noise4;
d4=abs(y4-xcomp);
figure(3);
hist(d3,1000);
grid on;
figure(4);
hist(d4,1000);
grid on;

% sum(d4.^2)/length(d4)
