clc;
clear;

isi=[0.28,1,0.28,4.3,7];
order=7;
snr=10;
len=1000;

M=(order-1)/2;
N=len+length(isi)-1;
e=zeros(1,len);
error=e;
out=zeros(1,N);
number=0;
r=zeros(order,order);

irt=1;
for j=1:irt
    theta=(rand-0.5)*2*pi;
    A=10;
    t=1:len;
    fc=1000;
    x_re=A*cos(2*pi*fc*t+theta);
    theta_im=(rand-0.5)*2*pi;
    A_im=5;
    fc_im=500;
    x_im=A_im*cos(2*pi*fc_im*t+theta_im);
    x=x_re.*exp(i*x_im)
    for p=1:order
        for q=p:order
            r(p,q)=A^2/2*cos(2*pi*fc*abs(p-q));
            r(q,p)=r(p,q);
        end
    end
    [V,D]=eig(r);
    lmd_max=max(max(D));
    mu=1/lmd_max;
        
    d1=zeros(1,len);
	noise=rand(1,N)/10.^(snr/10);
	y=conv(isi,x)+noise;
	wk=zeros(1,order);
	for n=1:len-order+1
		x1=x(n:n+order-1);
		d1(n)=wk*x1';
		e(n)=y(n)-d1(n);  
		wk=wk+mu*e(n)*x1;
    end
end
plot(abs(e))


