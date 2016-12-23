clc;
clear;

isi=[0.28,1,0.28];
order=5;
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
for i=1:irt
    theta=(rand-0.5)*2*pi;
    A=10;
    t=1:len;
    x=A*cos(0.1*t+theta); 
%     pn=rand(1,len);
%     sys=tf(1,[2 0 1]);
%     x=(lsim(sys,pn,1:len))';
    for p=1:order
        for q=p:order
            r(p,q)=A^2/2*cos(0.1*abs(p-q));
            r(q,p)=r(p,q);
        end
    end
    
    [V,D]=eig(r);
    lmd_max=max(max(D));
    mu=1/lmd_max;
    %mu=0.001;
         
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
plot(e)

% sys1=tf(isi,[1 0 0]);
% sys2=tf(wk,[1 0 0 0 0]);

% pause
% figure(1)
% freqz(isi,[1 0 0]);
% figure(2)
% freqz(fliplr(wk),[1 0 0 0 0]);

%%%%%% jianyan
% pause
% theta=(rand-0.5)*2*pi;
% x=A*cos(0.1*t+theta); 
% y=conv(isi,x)+noise;
% d1=zeros(1,length(y));
% for n=1:len-order+1
%     x1=x(n:n+order-1);
%     d1(n)=wk*x1';
%     e(n)=y(n)-d1(n);  
% end
% plot(1:1000,e(1:1000),'r',1:1000,y(1:1000),'g',1:1000,d1(1:1000),'b')
% title('ºÏ—È') 



