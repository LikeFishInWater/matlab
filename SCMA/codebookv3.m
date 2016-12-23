%reference£º Multi-dimensional SCMA Codebook Design Based on Constellation Rotation and Interleaving
% M:number of constellation points
% N:number of demension

%% mother constellation
M=4;
N=2;
S0=zeros(N,M);
for m=1:M
    S0(1,m)=(2*m-1-M)*(1+1i);
end
for n=2:N
    theta=(n-1)*pi/M/N;
    S0(n,:)=theta*S0(1,:);
end
% [temp,alpha]=sort(rand(1,M));
alpha=[3 1 4 2];
MC=S0;
for n=2:2:N
    MC(n,:)=S0(n,alpha);
end
%% rotation
eu=randi(10);
df=3;
J=6;
K=4;
for u=1:df
    fai(u)=(u-1)*2*pi/M/df+eu*2*pi/M;
end
F=[1 1 0 0;
    1 0 1 0;
    1 0 0 1 ;
    0 1 0 1 ;
    0 1 1 0;
    0 0 1 1;]';
F0=[fai(1) fai(2) 0 0;
    fai(2) 0 fai(3) 0;
    fai(3) 0 0 fai(1);
    0 fai(1) fai(2) 0;
    0 fai(3) 0 fai(2);
    0 0 fai(1) fai(3);]';
V(:,:,1)=[1 0;0 1;0 0;0 0];
V(:,:,2)=[1 0;0 0;0 1;0 0];
V(:,:,3)=[1 0;0 0;0 0;0 1];
V(:,:,4)=[0 0;1 0;0 0;0 1];
V(:,:,5)=[0 0;1 0;0 1;0 0];
V(:,:,6)=[0 0;0 0;1 0;0 1];
S=zeros(K,J,M);
for j=1:J
    for m=1:M
        S(:,j,m)=V(:,:,j)*MC(:,m).*exp(1i*F0(:,j));
    end
end
S=S/sqrt(17.3);
save Sv3.mat S F;

