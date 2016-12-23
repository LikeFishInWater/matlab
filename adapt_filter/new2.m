%y(n)=sum w(k)*n(n-k); k=0 1 2 â€¦â? 63
%e(n)=d(n)-y(n);
%w(k+1)=w(k)+2*mu*e(n)*
clc;
clear;
isi=[0.28,1,0.28];
order=63;
snr=30;
len=1000;
mu=0.02;

M=(order-1)/2;
N=len+length(isi)-1;
e=zeros(1,N);
error=e;
out=zeros(1,N);
number=0;
for i=1:1000
	x=sign(rand(1,len)-0.5);
	noise=rand(1,N)/10.^(snr/10);
	y=conv(isi,x)+noise;
	wk=zeros(1,order);
	for n=order:N-M+1
		y1=y(n+M-1:-1:n-M-1);
		d1=wk*y1';
		e(n)=x(n-2)-d1;
		wk=wk+mu*e(n)*y1;
		e(n)=10*log10(abs(e(n)));
	end
	error=error+e;
end
error=error(order:N-M+1)/100;

%%% help lsim