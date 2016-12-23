clc;
clear;

%%%%%%%%%%%
isi=[0.28,1,0.28];
snr=10;
len=1000; %—µ¡∑–Ú¡–
order=4;
M=order;
N=2*order+2;

delta=0.02;
lmd=0.98;

%%%%%%%%%%%%
theta=(rand-0.5)*2*pi;
A=10;
t=1:len;
x=A*cos(0.1*t+theta);

d1=zeros(1,len);
noise=rand(1,len+length(isi)-1)/10.^(snr/10);
y=conv(isi,x)+noise;

%%%%%%%%%%% rls¬À≤®
Wk=zeros(M,1);
P=delta*eye(M,M);
for n=1:len-M+1
    U=x(n:n+M-1)';
    d=y(n+M-1);
    Pai=P*U;
    K=Pai/(lmd+U'*Pai);
    ksi(n)=d-Wk'*U;
    Wk=Wk+K*conj(ksi(n))';
    P=1/lmd*P-1/lmd*K*U'*P;   
end
plot(ksi)


%%%%%%%%%%%  ºÏ—È
% pause
% B=20;
% theta=(rand-0.5)*2*pi;
% x=B*cos(0.3*t+theta); 
% y=conv(isi,x)+noise;
% d1=zeros(1,length(y));
% Wk=zeros(M,1);
% P=delta*eye(M,M);
% for n=1:len-M+1
%     U=x(n:n+M-1)';
%     d1(n+len-M)=Wk'*U;
% end
% e(1:1000)=y(1:1000)-d1(1:1000);
% plot(ksi)
% plot(1:1000,e(1:1000),'r',1:1000,y(1:1000),'g',1:1000,d1(1:1000),'b')



