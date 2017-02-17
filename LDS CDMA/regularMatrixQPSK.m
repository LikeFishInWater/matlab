clear;
%% qpsk regular-matrix lds-cdma uncoded
EsNo=11:20;
ite=7;
%% read F matrix
BinaryFile = 'MatszhN16_12_3_4_sample.txt';
fid=fopen(BinaryFile);
line=str2num(fgetl(fid));
K=line(1);
N=line(2); 
line=str2num(fgetl(fid));
dv=line(1);
dc=line(2);
fgetl(fid);fgetl(fid);
H_row_id = [];
H_col_id = [];
F=zeros(N,K);
for col = 1:K
    row_id = str2num(fgetl(fid));
    F(row_id,col)=1;
end
fclose(fid);
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
    for n=1:K
        if F(m,n)==1
            S(m,n)=sqrt(1/dc)*exp(i*pi/K/2*n);
        end
    end
end
H=S;

x_decode=zeros(K,1);
decode=zeros(K,1);

hModulator = comm.QPSKModulator('BitInput',true);
hModulator.PhaseOffset = 0;
status=zeros(dc,4^dc);
for ii=1:4^dc
    temp=ii-1;
    for jj=1:(dc-1)
        status(jj,ii)=mod(temp,4);
        temp=(temp-status(jj,ii))/4;
    end
    status(dc,ii)=temp;
%     status(1,ii)=mod(ii-1,4);
%     temp=((ii-1)-status(1,ii))/4;
%     status(2,ii)=mod(temp,4);
%     temp=(temp-status(2,ii))/4;
%     status(3,ii)=mod(temp,4);
%     temp=(temp-status(3,ii))/4;
%     status(4,ii)=(temp-status(3,ii))/4;
end
status=status+1;
status0=zeros(dc,4^dc);
for ii=1:4^dc
    for j=1:dc
        status0(j,ii)=exp((status(j,ii))*pi/2*1i);
    end
end
tic;


for ir=1:length(EsNo)
    snr=EsNo(ir);
    flag=1;
    fe(ir)=0;
    be(ir)=0;
    total(ir)=0;
    user_error=zeros(1,K);
    while flag
        total(ir)=total(ir)+1;
        x=randi([0 1],2*K,1);  
        modData = step(hModulator, x);
        temp=H*modData; 
        y=awgn(temp,snr);
        sigma=sqrt(1/2*1/(10^(snr/10)));
        %%%%%%%%%%%%%%% MPA decoder %%%%%%%%%%%%%%%%%%%%%
        D=zeros(N,K,4);
        for m=1:K
            for nn=1:dv
                D(v(nn,m),m,:)=0.25;
            end
        end
        for j=1:ite
            %%%%%%%%%%%%%%%%%Uk to Cn
            Dtemp=D;
            for m=1:K
                for nn=1:dv
                    Dtemp(v(nn,m),m,:)=0.25;
                    for tt=1:dv
                        if tt~=nn
                            Dtemp(v(nn,m),m,:)=Dtemp(v(nn,m),m,:).*D(v(tt,m),m,:);
                        end
                    end
                    Dtemp(v(nn,m),m,:)=Dtemp(v(nn,m),m,:)./sum(Dtemp(v(nn,m),m,:));
                end
            end
            D=Dtemp;    
            %%%%%%%%%%%%%Cn to Uk     
            Dtemp=zeros(N,K,4);
            for m=1:N
                for n=1:dc
                    for t0=1:4^dc
                        temp1=1/D(m,h(m,n),status(n,t0));
                        temp2=y(m);
                        for tt=1:dc
                            temp1=temp1*D(m,h(m,tt),status(tt,t0));
                            temp2=temp2-status0(tt,t0)*H(m,h(m,tt));
                        end
                        Dtemp(m,h(m,n),status(n,t0))= Dtemp(m,h(m,n),status(n,t0))+temp1*(1/(sqrt(2*pi)*sigma)*exp(-1/2/sigma^2*(abs(temp2))^2));
%                         Dtemp(m,h(m,n),status(n,t0))= Dtemp(m,h(m,n),status(n,t0))+...
%                             D(m,h(m,1),status(1,t0))*D(m,h(m,2),status(2,t0))*D(m,h(m,3),status(3,t0))*D(m,h(m,4),status(4,t0))/D(m,h(m,n),status(n,t0))...
%                             *(1/(sqrt(2*pi)*sigma)...
%                             *exp(-1/2/sigma^2*(abs(y(m)-status0(1,t0)*H(m,h(m,1))-status0(2,t0)*H(m,h(m,2))-status0(3,t0)*H(m,h(m,3))-status0(4,t0)*H(m,h(m,4))))^2));
                    end         
                end
            end
            D=Dtemp;
            %%%%%%%%%%%%%µü´úÖÐ¼ì²â  
        end
        for m=1:K
            temp=D(v(1,m),m,:);
            for nn=2:dv
                temp=temp.*D(v(nn,m),m,:);
            end
            mult(m,:)=temp;
            [~,decode(m)]=max(mult(m,:));
            switch decode(m)
                case 1,decodeBit([2*m-1,2*m])=[0 1];
                case 2,decodeBit([2*m-1,2*m])=[1 1];
                case 3,decodeBit([2*m-1,2*m])=[1 0];
                case 4,decodeBit([2*m-1,2*m])=[0 0];
            end
        end
        for p=1:2*K
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
    ratio_be(ir)=be(ir)/total(ir)/K
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
