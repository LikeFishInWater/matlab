%% lms algorithm
clear;
scale=2048;
%% noise
sigma_2=1;
sigma=sqrt(sigma_2);
noise=sigma*randn(scale,1);
%% e(n) to x(n),H(z)
a1=-1.6;
a2=0.8;
x(1)=noise(1);
x(2)=noise(2)-a1*x(1);
for n=3:scale
    x(n)=noise(n)-a1*x(n-1)-a2*x(n-2);
end
%% AF lms
mu=0.02;
order=2;
wk=zeros(scale-order,order);
for n=1:scale-order
    x1=x(n:n+order-1);
    d1(n)=wk(n,:)*x1';
    e(n+order)=x(n+order)-d1(n);  
    wk(n+1,:)=wk(n,:)+mu*e(n+order)*x1;
end
for n=1:2045
    AE(n)=sum(e(n:n+3))/4;
end

figure(1);
hold on;
title('随迭代次数变换的性能曲线')
plot(abs(e'));

grid on;
hold off;
figure(2);
hold on;
title('系数收敛曲线')
plot(wk(:,1),'r-');
plot(wk(:,2),'g-');
legend('a2收敛曲线','a1收敛曲线');

grid on;
hold off;
