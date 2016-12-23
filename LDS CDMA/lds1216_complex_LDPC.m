clear;
%% 
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
            S(m,n)=sqrt(1/4)*exp(i*pi/16*n);
        end
        % S(m,h(m,n))=sqrt(1/4)*exp(i*pi/2*n);
    end
end
H=S;
status=zeros(4,16);
for t=1:16
    a=de2bi(t-1,4);
    a=a*2-1;
    status(:,t)=a;
end

%% ldpc module 
dataLen=2000;
LdpcLen=4000;
BinaryFile = 'MN_4000.txt';
fid=fopen(BinaryFile);
line=str2num(fgetl(fid));
Nn=line(1);
Mm=line(2); 
Rr = (Nn - Mm) / Nn;
fgetl(fid);fgetl(fid);fgetl(fid);
H_row_id = [];
H_col_id = [];
for col = 1:Nn
    row_id = str2num(fgetl(fid));
    row_id(row_id == 0) = [];
    H_row_id = [H_row_id, row_id];
    col_id = col * ones(1,length(row_id));
    H_col_id=[H_col_id, col_id];
end
fclose(fid);
Hldpc = sparse(H_row_id, H_col_id, ones(1,length(H_row_id)));
Hldpc = (Hldpc >= 1);
clear H_row_id;clear H_col_id;
hEnc = comm.LDPCEncoder('ParityCheckMatrix' ,Hldpc);
hDec = comm.LDPCDecoder('ParityCheckMatrix' , Hldpc,'OutputValue','Whole codeword','DecisionMethod','Hard decision','IterationTerminationCondition','Parity check satisfied','FinalParityChecksOutputPort',true,'MaximumIterationCount',50);
%% 
LdsDecode=zeros(LdpcLen,16);
LdsDecLLR=zeros(LdpcLen,16);
tic;
EsNo=0.5:0.2:1.5;
for ir=1:length(EsNo)
    snr=EsNo(ir)
    sigma=sqrt(1/2*1/(10^(snr/10)));
    flag=1;
    errBit=0;
    errFrm=0;
    LDSerrBit=0;
    total=0;
    while flag
        total=total+1;
        for ii=1:16
            x(:,ii)=randi([0 1],dataLen,1);
            xLdpc(:,ii) = step(hEnc,x(:,ii));
        end
        for ii=1:LdpcLen
            temp=H*(2*xLdpc(ii,:)'-1); 
            y=awgn(temp,snr);
            %%%%%%%%%%%%%%% MPA decoder %%%%%%%%%%%%%%%%%%%%%
            D=zeros(12,16);
            for j=1:7
                %%%%%%%%%%%%%%%%%Uk to Cn
                for m=1:16;
                    he=0;
                    for n=1:3
                        he=he+D(v(n,m),m);
                    end
                    for n=1:3
                        D(v(n,m),m)=he-D(v(n,m),m);
                    end
                end        
                %%%%%%%%%%%%%Cn to Uk                 
                D_temp=D;
                for n=1:12  
                    for c=1:4
                        rpc=1;rnc=1;
                        for t=1:16
                            if(status(c,t)==1)
                                rp(rpc)=status(1,t)/2*D(n,h(n,1))+status(2,t)/2*D(n,h(n,2))+status(3,t)/2*D(n,h(n,3))+status(4,t)/2*D(n,h(n,4))-status(c,t)/2*D(n,h(n,c))-1/(2*sigma^2)*(abs(y(n)-status(1,t)*H(n,h(n,1))-status(2,t)*H(n,h(n,2))-status(3,t)*H(n,h(n,3))-status(4,t)*H(n,h(n,4))))^2;
                                rpc=rpc+1;
                            else
                                rn(rnc)=status(1,t)/2*D(n,h(n,1))+status(2,t)/2*D(n,h(n,2))+status(3,t)/2*D(n,h(n,3))+status(4,t)/2*D(n,h(n,4))-status(c,t)/2*D(n,h(n,c))-1/(2*sigma^2)*(abs(y(n)-status(1,t)*H(n,h(n,1))-status(2,t)*H(n,h(n,2))-status(3,t)*H(n,h(n,3))-status(4,t)*H(n,h(n,4))))^2;
                                rnc=rnc+1;
                            end
                        end
                        % D_temp(n,h(n,c))=log(sum(exp(rp))/sum(exp(rn)));
                        ptem=rp(1);
                        ntem=rn(1);
                        for t=2:8
                            ptem=max(ptem,rp(t));
                            ntem=max(ntem,rn(t));
                        end
                        D_temp(n,h(n,c))=ptem-ntem;
                    end
                end
                D=D_temp;
            end
            for k=1:16
                he=0;
                for n=1:12
                    he=he+D(n,k);
                end
                LdsDecLLR(ii,k)=he;
                if he>0
                    LdsDecode(ii,k)=1;
                else
                    LdsDecode(ii,k)=0;
                end
            end          
        end
        %%%%%%%%%ldpc decode
        for ii=1:16
            LdsDec(:,ii)=step(hDec,-LdsDecLLR(:,ii));
        end
        for ii=1:16
            for jj=1:LdpcLen
                if LdsDecode(jj,ii)==xLdpc(jj,ii)
                    LDSerrBit=LDSerrBit;
                else
                    LDSerrBit=LDSerrBit+1;
                end
            end
        end
        for ii=1:16
            for jj=1:dataLen
                if(LdsDec(jj,ii)==x(jj,ii))
                    errBit=errBit;
                else
                    errBit=errBit+1;
                end
            end
            if(LdsDec(1:dataLen,ii)==x(:,ii))
                errFrm=errFrm;
            else 
                errFrm=errFrm+1;
            end
        end
        if errFrm>=20
            flag=0;
        end
        %%%%%%%%%ldpc decode end
    end
    ber(ir)=errBit/total/16/dataLen
    fer(ir)=errFrm/total/16
    LDSber(ir)=LDSerrBit/total/16/LdpcLen
toc;
end
% semilogy(EsNo,ber,'bo-');
% legend('lds1216');hold on;
% xlabel('snr');ylabel('ber');title('lds-cdma12*16');
% axis([(1+ir_snr) (ir_max+ir_snr) 1e-5 1e-1]);
% grid on;
% hold off;
% str=['figure',num2str(sav)];
% str2=['S_mat',num2str(sav),'.mat'];
% print(1,'-dpng',str);
% save(str2,'S');
% end
% ber(sav)=ratio_be;
% end
% save lds1216_complex_LDPC.mat;


