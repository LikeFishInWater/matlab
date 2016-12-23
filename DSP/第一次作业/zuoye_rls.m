%% rls algorithm
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
%% AF rls
delta=0.02;
lmd=0.9;
order=2;
Wk=zeros(order,scale-order);
P=1/delta*eye(order,order);
for n=1:scale-order
    U=((fliplr(x(n:n+order-1)))');
    d=x(n+order);
    Pai=P*U;
    K=Pai/(lmd+U'*Pai);
    ksi(n+order)=d-(Wk(:,n)'*U);
    Wk(:,n+1)=Wk(:,n)+K*(ksi(n+order));
    P=1/lmd*P-1/lmd*K*U'*P;   
end
plot(abs(ksi));
axis([0 200 0 10])
figure(1);
hold on;
title('随迭代次数变换的性能曲线')
plot(abs(ksi));
axis([0 2048 0 10])
grid on;
hold off;

figure(2);
hold on;
title('系数收敛曲线')
plot(Wk(1,:),'r-');
plot(Wk(2,:),'g-');
legend('a2收敛曲线','a1收敛曲线');
grid on;
hold off;