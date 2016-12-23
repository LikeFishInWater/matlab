clc;
clear;

%%%%%%%%%%%
isi=[0.28,1,0.28];
snr=10;
len=200; %—µ¡∑–Ú¡–
order=7;
M=order;
N=2*order+2;

delta=0.1;
lmd=0.9;

%%%%%%%%%%%%
theta=(rand-0.5)*2*pi;
A=1;
t=1:len;
x=A*(cos(0.1*t+theta)+ 1i*sin(0.1*t+theta));%

d1=zeros(1,len);
noisere=1/2*rand(1,len+length(isi)-1)/10.^(snr/10);
noiseim=1/2*rand(1,len+length(isi)-1)/10.^(snr/10);
noise=noisere+1i*noiseim;
y=conv(isi,x);

%%%%%%%%%%% rls¬À≤®
Wk=zeros(M,1);
P=1/delta*eye(M,M);
for n=1:len-M+1
    U=conj((fliplr(x(n:n+M-1)))');
    d=y(n+M-1);
    Pai=P*U;
    K=Pai/(lmd+U'*Pai);
    ksi(n+M-1)=d-(Wk'*U);
    Wk=Wk+K*conj(ksi(n+M-1));
    P=1/lmd*P-1/lmd*K*U'*P;   
end
plot(abs(ksi))




