clear;
load test.mat;
scale=length(x);

%% r(1)代表r0,依次类推
r=zeros(1,2*scale-1);
r_i=-scale+1:1:scale-1;
for ii=(scale-1):-1:1
    for jj=scale:-1:(1+ii)
        r(scale-ii)=r(scale-ii)+x(jj)'*x(jj-ii);
    end
    r(scale-ii)=r(scale-ii)/(scale);
end
for ii=1:scale
    for jj=1:(scale-ii+1)
        r(scale+ii-1)=r(scale+ii-1)+x(jj)'*x(jj+ii-1);
    end
    r(scale+ii-1)=r(scale+ii-1)/(scale);
end
wb=-pi:0.001:pi;
fb=wb/2/pi;
%% 理论值
% figure(1)
% sigma_2=0.01;
% sigma=sqrt(sigma_2);
% u=sigma*randn(10000,1);
% v=sigma*randn(10000,1);
% y=u+1j*v;
% dB=[64 54 2 30];
% A=sqrt(10.^(dB/10)*0.02);
% f=[0.15 0.16 0.252 -0.16];
% for ii=1:10000
%     x0(ii)=y(ii)+sum(A.*exp(1j*2*pi*f*ii));
% end
% dtft_x=zeros(1,length(wb));
% for ii=1:10000
%     dtft_x=dtft_x+x0(ii)*exp(-1j*ii.*wb);
% end
% pb=1/10000*(abs(dtft_x).^2);
% pb=pb/sum(pb);
% semilogy(fb,pb);
% grid on;

%% 直接法（周期图法）：先对x求fft再平方
figure(2)
dtft_x=zeros(1,length(wb));
for ii=1:scale
    dtft_x=dtft_x+x(ii)*exp(-1j*ii.*wb);
end
pb=1/scale*(abs(dtft_x).^2);
pb=pb/max(pb);
pb_db=10*log10(pb);
plot(fb,pb_db);
%semilogy(fb,pb);
xlabel('w/(2pi)');ylabel('归一化P(w)/dB');
title('直接法（周期图法）');
axis([-0.5 0.5 -50 0]);
grid on;
%% 平均周期图，四段，每段32，无重叠，汉明窗
figure(3);
num=4;
subScale=scale/num;
w=hamming(subScale);
pc=zeros(num,length(wb));
dtft_x=zeros(num,length(wb));
for ii=1:num
    for jj=1:subScale
        dtft_x(ii,:)=dtft_x(ii,:)+x((ii-1)*subScale+jj)*exp(-1j*(jj)*wb);
    end
    pc(ii,:)=1/subScale*(abs(dtft_x(ii,:)).^2);
    pc(ii,:)=pc(ii,:)/max(pc(ii,:));
end
pc_average=sum(pc,1)/num;
pc_average_db=10*log10(pc_average);
plot(fb,pc_average_db);
% semilogy(fb,pc_average);
xlabel('w/(2pi)');ylabel('归一化P(w)/dB');
title('平均周期图，四段，每段32，无重叠，汉明窗');
axis([-0.5 0.5 -50 0]);
grid on;
%% 平均周期图，每段32点，叠合16点，汉明窗
figure(4);
subScale=32;
overlap=16;
num=7;
pd=zeros(num,length(wb));
dtft_x=zeros(num,length(wb));
for ii=1:num
    for jj=1:subScale
        dtft_x(ii,:)=dtft_x(ii,:)+x((ii-1)*(subScale-overlap)+jj)*exp(-1j*jj*wb);
    end
    pd(ii,:)=1/subScale*(abs(dtft_x(ii,:)).^2);
    pd(ii,:)=pd(ii,:)/max(pd(ii,:));
end
pd_average=sum(pd,1)/num;
pd_average_db=10*log10(pd_average);
plot(fb,pd_average_db);
% semilogy(fb,pd_average);
xlabel('w/(2pi)');ylabel('归一化P(w)/dB');
title('平均周期图，每段32点，叠合16点，汉明窗');
axis([-0.5 0.5 -50 0]);
grid on;
%% 间接法：先求自相关函数再dtft,M=32
figure(5);
M=32;
r0=r((scale-M):(scale+M));
dtft_r0=zeros(1,length(wb));
for ii=1:(2*M+1)
    dtft_r0=dtft_r0+r0(ii)*exp(-1j*(ii-M-1)*wb);
end
pe=abs(dtft_r0);
pe=pe/max(pe);
pe_db=10*log10(pe);
plot(fb,pe_db)
% semilogy(fb,abs(pe));
xlabel('w/(2pi)');ylabel('归一化P(w)/dB');
title('间接法，M=32');
axis([-0.5 0.5 -50 0]);
grid on;
%% 间接法：先求自相关函数再dtft,M=16，汉明窗
figure(6);
M=16;
w=hamming(2*M+1);
re=r((scale-M):(scale+M)).*w';
dtft_r0=zeros(1,length(wb));
for ii=1:(2*M+1)
    dtft_r0=dtft_r0+re(ii)*exp(-1j*(ii-M-1)*wb);
end
pf=abs(dtft_r0);
pf=pf/max(pf);
pf_db=10*log10(pf);
plot(fb,pf_db);
% semilogy(fb,pf);
xlabel('w/(2pi)');ylabel('归一化P(w)/dB');
title('间接法，M=16，汉明窗');
axis([-0.5 0.5 -50 0]);
grid on;