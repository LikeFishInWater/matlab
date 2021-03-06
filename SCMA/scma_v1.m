% clc;
clear;

M=4;
N=2;
K=4;
J=6;

F=[0 1 0 1;
    1 0 1 0;
    1 1 0 0;
    0 0 1 1;
    1 0 0 1;
    0 1 1 0;]';
S=zeros(K,J,M);
S(:,:,1)=[0,-0.1815-0.1318*i,0,0.7851;0.7851,0,-0.7815-0.1318*i,0;-0.6351+0.4615*i,0.1392-0.1759*i,0,0;
    0,0,0.7851,-0.0055-0.2242*i;-0.0055-0.2242*i,0,0,-0.6351+0.4615*i;0,0.7851,0.1392-0.1759*i,0;].';
S(:,:,2)=[0,-0.6351-0.4615*i,0,-0.2243;-0.2243,0,-0.6351-0.4615*i,0;0.1815-0.1318*i,0.4873-0.6156*i,0,0;
    0,0,-0.2243,-0.0193-0.7848*i;-0.0193-0.7848*i,0,0,0.1815-0.1318*i;0,-0.2243,0.4873-0.6156*i,0;].';
S(:,:,3)=[0,0.6351+0.4615*i,0,0.2243;0.2242,0,0.6351+0.4615*i,0;-0.1815+0.1318*i,-0.4873+0.6156*i,0,0;
    0,0,0.2243,0.0193+0.7848*i;0.0193+0.7848*i,0,0,-0.1815+0.1318*i;0,0.2243,-0.4873+0.6156*i,0;].';
S(:,:,4)=[0,0.1815+0.1318*i,0,-0.7851;-0.7851,0,0.1815+0.1318*i,0;0.6351-0.4615*i,-0.1392+0.1759*i,0,0;
    0,0,-0.7851,0.0055+0.2242*i;0.0055+0.2242*i,0,0,0.6351-0.4615*i;0,-0.7851,-0.1392+0.1759*i,0;].';

% F=[1 1 0 0;1 0 1 0;1 0 0 1;0 1 1 0;0 1 0 1;0 0 1 1;]';
% R1=1;
% theta1=pi/3;
% theta2=pi/3*2;
% theta3=pi;
% alpha=3;
% beta=1/0.63;
% R2=R1*beta;
% MC=[alpha*R1,R1,-R1,-alpha*R1;
%     -R2,alpha*R2,-alpha*R2,R2];
% delta1=[exp(1i*theta1),0;0,exp(1i*theta2)];
% delta2=[1,0;0,1];
% delta3=[0,1;1,0]*[exp(1i*theta1),0;0,exp(1i*theta3)];
% delta4=[1,0;0,exp(1i*theta2)];
% delta5=[0,1;1,0];
% delta6=[1,0;0,exp(1i*theta3)];
% V1=[1 0;0 1;0 0;0 0];
% V2=[1 0;0 0;0 1;0 0];
% V3=[1 0;0 0;0 0;0 1];
% V4=[0 0;1 0;0 1;0 0];
% V5=[0 0;1 0;0 0;0 1];
% V6=[0 0;0 0;1 0;0 1];
% C1=V1*delta1*MC;
% C2=V2*delta2*MC;
% C3=V3*delta3*MC;
% C4=V4*delta4*MC;
% C5=V5*delta5*MC;
% C6=V6*delta6*MC;
% S=zeros(K,J,M);
% S(:,1,:)=C1;
% S(:,2,:)=C2;
% S(:,3,:)=C3;
% S(:,4,:)=C4;
% S(:,5,:)=C5;
% S(:,6,:)=C6;
% S=S/sqrt(26.32);
h=zeros(4,3);
for m=1:4        
    h(m,:)=find(F(m,:));
end
v=zeros(2,6);
for m=1:6
	v(:,m)=find(F(:,m));
end
pp=zeros(64800/2,M);
for jj=1:64800/2
    for mm=1:M
        pp(jj,mm)=1/M;
    end
end
status=zeros(2,16);
for n=1:16
    status(1,n)=mod(n-1,4);
    status(2,n)=(n-1-status(1,n))/4;
end
status=status+1;

EsNo=10:20;
for mm=1:length(EsNo)
    flag=1;
    snr=EsNo(mm);
    fe(mm)=0;
    be(mm)=0;
    bbe=0;
    total(mm)=0;
    while flag
        total(mm)=total(mm)+1;
        for ii=1:J
            u(ii)=randi([1,M]);
        end
        % u=[ 1     4     2     2     2     3];
        x=zeros(K,1);
        for ii=1:J
            x=x+S(:,ii,u(ii));
        end
        % noise=zeros(K,1);
        % y=x+noise;
        y=awgn(x,snr);
        sigma=sqrt(1/2*10^(-snr/10));
        %%%%%%%%% decoding %%%%%%%%%%%
        D=zeros(K,J,M);
        for m=1:4
            for n=1:3
                D(m,h(m,n),:)=[1 1 1 1];
            end
        end
        T=7;
        for t=1:T
            Dtemp=D;
            for m=1:J
                for n=1:M
                    Dtemp(v(1,m),m,n)=Dtemp(v(1,m),m,n).*pp(((mm-1)*6+m),n);
                    Dtemp(v(2,m),m,n)=Dtemp(v(2,m),m,n).*pp(((mm-1)*6+m),n);
                end
            end
            for m=1:J
                D(v(1,m),m,:)=Dtemp(v(2,m),m,:)/sum(Dtemp(v(2,m),m,:));
                D(v(2,m),m,:)=Dtemp(v(1,m),m,:)/sum(Dtemp(v(1,m),m,:));
            end
            Dtemp=D;
            for k=1:4
                j=1;
                for m=1:4
                    for n=1:16	
                        p(n)=D(k,h(k,2),status(1,n))*D(k,h(k,3),status(2,n))*exp(-1/(2*sigma^2)*(abs(y(k)-S(k,h(k,1),m)-S(k,h(k,2),status(1,n))-S(k,h(k,3),status(2,n))))^2);
                        Dtemp(k,h(k,j),m)=sum(p);
                    end
                end
                j=2;
                for m=1:4
                    for n=1:16
                        p(n)=D(k,h(k,1),status(1,n))*D(k,h(k,3),status(2,n))*exp(-1/(2*sigma^2)*(abs(y(k)-S(k,h(k,2),m)-S(k,h(k,1),status(1,n))-S(k,h(k,3),status(2,n))))^2);
                        Dtemp(k,h(k,j),m)=sum(p);
                    end
                end
                j=3;
                for m=1:4
                    for n=1:16
                        p(n)=D(k,h(k,1),status(1,n))*D(k,h(k,2),status(2,n))*exp(-1/(2*sigma^2)*abs(y(k)-S(k,h(k,3),m)-S(k,h(k,1),status(1,n))-S(k,h(k,2),status(2,n)))^2);
                        Dtemp(k,h(k,j),m)=sum(p);
                    end
                end
            end
            D=Dtemp;
        end

        De=zeros(1,J,M);
        for j=1:J
            De(1,j,:)=D(v(1,j),j,:).*D(v(2,j),j,:);
            De(1,j,:)=De(1,j,:)/sum(De(1,j,:));
        end

        for j=1:J
            [b,decode(j)]=max(De(1,j,:));
        end
        for j=1:J
            if decode(j)==u(j)
                be(mm)=be(mm);
            else
                be(mm)=be(mm)+1;
            end
        end
        for j=1:J
            if abs(decode(j)-u(j))==2
                bbe=bbe+2;
            else
                if decode(j)==u(j)
                    bbe=bbe;
                else
                    bbe=bbe+1;
                end
            end
        end
        if decode==u
            fe(mm)=fe(mm);
        else
            fe(mm)=fe(mm)+1
        end
        if fe(mm)==30
            flag=0;
        end
    end
    ratio_fe(mm)=fe(mm)/total(mm)
    ratio_be(mm)=be(mm)/total(mm)/6
    ratio_bbe(mm)=bbe/total(mm)/12
end
semilogy(5:10,ratio_bbe);
grid on;
xlabel('snr(Es/N0)');ylabel('ber');
title('scma-0406-complex-uncoded');
axis([5 10 1e-2 1e0]);
save 
