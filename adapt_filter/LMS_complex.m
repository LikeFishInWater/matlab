clc;
clear;

isi=[0.28,1,0.28,4.3,7]; %信道参数
order=63; %阶数
snr=10; %信噪比
len=2000; %输入长度（与训练序列长度有关）
adjust=1; %步长调整系数
% accuracy=64; %计算精度
% 
% digits(accuracy);

N=len+length(isi)-1;
r=zeros(order,order);

%输入序列产生
theta=(rand-0.5)*2*pi;
t=1:len;
x=cos(t+theta) + 1i*sin(t+theta); %exp(i(t+theta))

%计算输入信号自相关矩阵r,获得步长mu
for p=1:order
    for q=1:order
        r(p,q)=exp(1i * (p-q));
    end
end
[V,D]=eig(r);
lmd_max=max(max(D));
mu=1/lmd_max * adjust;

%lms算法
d1=zeros(1,len);
noise=  0.1 * (randn(1,N)/ 10.^(snr/10) + 1i * randn(1,N)/ pi * 10.^(snr/10));
y=conv(isi,x)+noise;
wk=zeros(1,order);
for n=order:len  
    x1=x(n - order + 1: n);
    d1(n)=x1*wk';
    e(n)=y(n)-d1(n);  
    wk=wk+mu*conj(e(n))*x1;
end
sum = e(1000:2000)*e(1000:2000)';
n = noise(1000:2000) * noise(1000:2000)';
k = (sum - n)/1000;
plot(abs(e));
xlabel ('训练长度');ylabel ('误差');
axis([0 300 0 10]);
STR=sprintf('%s%f','步长系数为',adjust);
title (STR);

