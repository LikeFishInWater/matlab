clear;
%indicator matrix
F=[0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0;
    0,0,0,1,0,0,1,0,0,0,0,0,1,0,1,0;
    0,1,0,0,1,0,1,0,0,1,0,0,0,0,0,0;
    0,0,0,1,0,1,0,0,0,0,1,0,0,1,0,0;
    0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,1;
    1,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0;
    0,0,0,1,0,0,0,1,1,0,0,1,0,0,0,0;
    0,1,0,0,0,1,0,0,0,0,0,1,0,0,0,1;
    1,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1;
    0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0;
    0,1,0,0,0,0,0,0,1,0,1,0,1,0,0,0;
    1,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0;];
%spreading signature matrices
h=zeros(12,4);
for m=1:12          
    h(m,:)=find(F(m,:));
end
v=zeros(3,16);
for m=1:16
    v(:,m)=find(F(:,m));
end
S=zeros(12,16);
for m=1:12
    for n=1:16
        if F(m,n)==1
            S(m,n)=sqrt(1/4)*exp(i*pi/32*n);
        end
        % S(m,h(m,n))=sqrt(1/4)*exp(i*pi/2*n);
    end
end
H=S;

x_decode=zeros(16,1);
decode=zeros(16,1);

hModulator = comm.QPSKModulator('BitInput',true);
hModulator.PhaseOffset = 0;
status=zeros(4,256);
for ii=1:256
    status(1,ii)=mod(ii-1,4);
    temp=((ii-1)-status(1,ii))/4;
    status(2,ii)=mod(temp,4);
    temp=(temp-status(2,ii))/4;
    status(3,ii)=mod(temp,4);
    status(4,ii)=(temp-status(3,ii))/4;
end
status=status+1;
status0=zeros(4,256);
for ii=1:256
    for j=1:4
        status0(j,ii)=exp((status(j,ii))*pi/2*1i);
    end
end
tic;
EsNo=11:20;
for ir=1:length(EsNo)
    snr=EsNo(ir);
    flag=1;
    fe(ir)=0;
    be(ir)=0;
    total(ir)=0;
    user_error=zeros(1,16);
    while flag
        total(ir)=total(ir)+1;
        x=randi([0 1],32,1);  
        modData = step(hModulator, x);
        temp=H*modData; 
        y=awgn(temp,snr);
        sigma=sqrt(1/2*1/(10^(snr/10)));
        %%%%%%%%%%%%%%% MPA decoder %%%%%%%%%%%%%%%%%%%%%
        D=zeros(12,16,4);
        for m=1:16
            D(v(1,m),m,:)=0.25;
            D(v(2,m),m,:)=0.25;
            D(v(3,m),m,:)=0.25;
        end
        for j=1:7
            %%%%%%%%%%%%%%%%%Uk to Cn
            Dtemp=D;
            for m=1:16
                Dtemp(v(1,m),m,:)=D(v(2,m),m,:).*D(v(3,m),m,:)*0.25;
                Dtemp(v(1,m),m,:)=Dtemp(v(1,m),m,:)./sum(Dtemp(v(1,m),m,:));
                Dtemp(v(2,m),m,:)=D(v(1,m),m,:).*D(v(3,m),m,:)*0.25;
                Dtemp(v(2,m),m,:)=Dtemp(v(2,m),m,:)./sum(Dtemp(v(2,m),m,:));
                Dtemp(v(3,m),m,:)=D(v(1,m),m,:).*D(v(2,m),m,:)*0.25;
                Dtemp(v(3,m),m,:)=Dtemp(v(3,m),m,:)./sum(Dtemp(v(3,m),m,:));
            end
            D=Dtemp;    
            %%%%%%%%%%%%%Cn to Uk     
            Dtemp=zeros(12,16,4);
            for m=1:12
                for n=1:4
                    for t0=1:256
                        Dtemp(m,h(m,n),status(n,t0))= Dtemp(m,h(m,n),status(n,t0))+...
                            D(m,h(m,1),status(1,t0))*D(m,h(m,2),status(2,t0))*D(m,h(m,3),status(3,t0))*D(m,h(m,4),status(4,t0))/D(m,h(m,n),status(n,t0))...
                            *(1/(sqrt(2*pi)*sigma)...
                            *exp(-1/2/sigma^2*(abs(y(m)-status0(1,t0)*H(m,h(m,1))-status0(2,t0)*H(m,h(m,2))-status0(3,t0)*H(m,h(m,3))-status0(4,t0)*H(m,h(m,4))))^2));
                    end         
                end
%                 n=1;
%                 for t=1:4
%                     for t0=1:64
%                         ptemp(t0)=D(m,h(m,2),status(1,t0))*D(m,h(m,3),status(2,t0))*D(m,h(m,4),status(3,t0))*(1/(sqrt(2*pi)*sigma)*exp(-1/2/sigma^2*(abs(y(m)-exp(pi/2*t*1i)*H(m,h(m,1))-status0(1,t0)*H(m,h(m,2))-status0(2,t0)*H(m,h(m,3))-status0(3,t0)*H(m,h(m,4))))^2));
%                     end
%                     Dtemp(m,h(m,n),t)=sum(ptemp);
%                 end
%                 n=2;
%                 for t=1:4
%                     for t0=1:64
%                         ptemp(t0)=D(m,h(m,1),status(1,t0))*D(m,h(m,3),status(2,t0))*D(m,h(m,4),status(3,t0))*(1/(sqrt(2*pi)*sigma)*exp(-1/2/sigma^2*(abs(y(m)-exp(pi/2*t*1i)*H(m,h(m,2))-status0(1,t0)*H(m,h(m,1))-status0(2,t0)*H(m,h(m,3))-status0(3,t0)*H(m,h(m,4))))^2));
%                     end
%                     Dtemp(m,h(m,n),t)=sum(ptemp);
%                 end
%                 n=3;
%                 for t=1:4
%                     for t0=1:64
%                         ptemp(t0)=D(m,h(m,1),status(1,t0))*D(m,h(m,2),status(2,t0))*D(m,h(m,4),status(3,t0))*(1/(sqrt(2*pi)*sigma)*exp(-1/2/sigma^2*(abs(y(m)-exp(pi/2*t*1i)*H(m,h(m,3))-status0(1,t0)*H(m,h(m,1))-status0(2,t0)*H(m,h(m,2))-status0(3,t0)*H(m,h(m,4))))^2));
%                     end
%                     Dtemp(m,h(m,n),t)=sum(ptemp);
%                 end
%                 n=4;
%                 for t=1:4
%                     for t0=1:64
%                         ptemp(t0)=D(m,h(m,1),status(1,t0))*D(m,h(m,2),status(2,t0))*D(m,h(m,3),status(3,t0))*(1/(sqrt(2*pi)*sigma)*exp(-1/2/sigma^2*(abs(y(m)-exp(pi/2*t*1i)*H(m,h(m,4))-status0(1,t0)*H(m,h(m,1))-status0(2,t0)*H(m,h(m,2))-status0(3,t0)*H(m,h(m,3))))^2));
%                     end
%                     Dtemp(m,h(m,n),t)=sum(ptemp);
%                 end
            end
            D=Dtemp;
            %%%%%%%%%%%%%µü´úÖÐ¼ì²â  
        end
        for m=1:16
            mult(m,:)=D(v(1,m),m,:).*D(v(2,m),m,:).*D(v(3,m),m,:);
            [~,decode(m)]=max(mult(m,:));
            switch decode(m)
                case 1,decodeBit([2*m-1,2*m])=[0 1];
                case 2,decodeBit([2*m-1,2*m])=[1 1];
                case 3,decodeBit([2*m-1,2*m])=[1 0];
                case 4,decodeBit([2*m-1,2*m])=[0 0];
            end
        end
        for p=1:32
            if decodeBit(p)==x(p)
                be(ir)=be(ir);
            else
                be(ir)=be(ir)+1;
            end
        end
        if decodeBit==x'
            fe(ir)=fe(ir);
        else
            fe(ir)=fe(ir)+1     
        end 
        if(fe(ir)==30)
            flag=0;
        end
    end
    ratio_fe(ir)=fe(ir)/total(ir)
    ratio_be(ir)=be(ir)/total(ir)/16
    toc;
end
semilogy((1+ir_snr):(ir_max+ir_snr),ratio_be,'bo-');
legend('lds1216');hold on;
xlabel('snr');ylabel('ber');title('lds-cdma12*16');
axis([(1+ir_snr) (ir_max+ir_snr) 1e-5 1e-1]);
grid on;
% hold off;
% str=['figure',num2str(sav)];
% str2=['S_mat',num2str(sav),'.mat'];
% print(1,'-dpng',str);
% save(str2,'S');
% end
% ber(sav)=ratio_be;
% end
save lds1216_complex.mat;


