clear;

%% read F matrix
F=[ 1 1 0 0;
    1 0 1 0;
    1 0 0 1;
    0 1 1 0;
    0 1 0 1;   
    0 0 1 1;]';
dc=3;dv=2;
K=6;N=4;
%% 
h=zeros(N,dc);
for m=1:N         
    h(m,:)=find(F(m,:));
end
v=zeros(dv,K);
for m=1:K
    v(:,m)=find(F(:,m));
end
S=zeros(N,K);
for m=1:N
    for n=1:dc
        S(m,h(m,n))=sqrt(1/dc)*exp(i*pi/dc/2*(n));
    end
end
H=S;

% hModulator = comm.QPSKModulator('BitInput',true);
% hModulator.PhaseOffset = 0;
scale=100000;
for ii=1:scale
    x=2*randi([0 1],K,1)-1;
    % modData = step(hModulator, x);
    temp=H*x; 
    average(ii)=sum(abs(temp).^2)/N;
end
sum(average)/scale
