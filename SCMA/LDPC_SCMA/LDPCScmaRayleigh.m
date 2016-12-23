clear;
%%%%%%%%%% scma module %%%%%%%%%%%%%
M=4;
N=2;
K=4;
J=6;
% F=[1 1 0 0
%     0 0 1 1
%     1 0 1 0
%     0 1 0 1
%     1 0 0 1
%     0 1 1 0]';
% load('S.mat');
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
h=zeros(4,3);
for m=1:4        
    h(m,:)=find(F(m,:));
end
v=zeros(2,6);
for m=1:6
	v(:,m)=find(F(:,m));
end
EsNo=10;

%%%%%%%%%%% ldpc module 
dataLen=2000;%%%%%%%%%%%%%%%%%%%
LdpcLen=4000;
R=1/2;
% hEnc = comm.LDPCEncoder();
% hDec = comm.LDPCDecoder('OutputValue','Whole codeword','DecisionMethod','Soft decision');
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
H = sparse(H_row_id, H_col_id, ones(1,length(H_row_id)));
H = (H >= 1);
clear H_row_id;clear H_col_id;
hEnc = comm.LDPCEncoder('ParityCheckMatrix' ,H);
hDec = comm.LDPCDecoder('ParityCheckMatrix' , H,'OutputValue','Whole codeword','DecisionMethod','Soft decision','MaximumIterationCount',50);

ModBits = 2;
hMod = comm.RectangularQAMModulator('ModulationOrder',2^ModBits,'BitInput',true,'NormalizationMethod', 'Average power');

for j=1:J
    [temp, alpha(:,j)] = sort(rand(1,LdpcLen));  
end
%%%%%%%%%%%%%%%%%
tic;
data=zeros(dataLen,J);
encodedData=zeros(LdpcLen,J);
IencodedData=encodedData;
modSignal=zeros(dataLen,J);
uerrBit=zeros(length(EsNo),J);
receivedBits=data;
itera=8;
a=zeros(2,16);
for i=1:16
    a(1,i)=mod(i-1,4);
    a(2,i)=(i-1-a(1,i))/4;
end
a=a+1;
Dinit=zeros(K,J,M);
for m=1:4
    for n=1:3
        Dinit(m,h(m,n),:)=1/4*[1 1 1 1];
    end
end
De=zeros(1,J,M); 
p=zeros(1,16);
yt=zeros(K,LdpcLen/2);
Isoft_demod_output=zeros(LdpcLen,J);
soft_demod_output=zeros(LdpcLen,J);
LdpcDecode=zeros(LdpcLen,J);
ILdpcDecode=zeros(LdpcLen,J);
outIteMax=2;
T=8;
for ii=1:length(EsNo)
    snr=EsNo(ii)
    sigma=sqrt(1/2*10^(-snr/10));
    flag=1;
    total=0;
    errBlc=0;
    errBit=0;
    %while flag 
    while flag
        pp=1/M*(zeros(LdpcLen/2,J,M)+1);  
        total=total+1;
        %%%%%%%%%% LDPC coding and modulation module %%%%%%%%%%%
        for j=1:J
            data(:,j) = randi([0 1], dataLen, 1);
            encodedData(:,j) = step(hEnc,data(:,j));
            IencodedData(alpha(:,j),j)=encodedData(:,j);
            for p1=1:LdpcLen/2
                modSignal(p1,j)=bin2dec(sprintf('%d',IencodedData(p1*2-1,j),IencodedData(p1*2,j)))+1;
            end
        end       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% transmition and receiving %%%%%%%%%%%
        for p1=1:LdpcLen/2
            x=zeros(K,1);
            for q1=1:J
                x=x+S(:,q1,modSignal(p1,q1));
            end
			ht = 1/sqrt(pi/2)*(randn(K,1) + 1i * randn(K,1));
			noise=sigma*(randn(K,1)+1i*randn(K,1));
            yt(:,p1)=x.*ht+noise;
        end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        toc
        for outIte=1:outIteMax
            %%%%%%%%%% SCMA coding and decoding module %%%%%%%%%%%%
            for p1=1:LdpcLen/2
                y=yt(:,p1);
                D=Dinit;
                for t=1:T
                    Dtemp=D;
                    for m=1:J
                        Dtemp(v(1,m),m,:)=Dtemp(v(1,m),m,:).*pp(p1,m,:);
                        Dtemp(v(2,m),m,:)=Dtemp(v(2,m),m,:).*pp(p1,m,:);
                    end
                    for m=1:J
                        D(v(1,m),m,:)=Dtemp(v(2,m),m,:)/sum(Dtemp(v(2,m),m,:));
                        D(v(2,m),m,:)=Dtemp(v(1,m),m,:)/sum(Dtemp(v(1,m),m,:));
                    end
                    Dtemp=D;
                    for k=1:4
                        for m=1:4
                            for n=1:16
                                p(n)=D(k,h(k,2),a(1,n))*D(k,h(k,3),a(2,n))*exp(-1/(2*sigma^2)*(abs(y(k)-S(k,h(k,1),m)-S(k,h(k,2),a(1,n))-S(k,h(k,3),a(2,n))))^2);
                            end
                            Dtemp(k,h(k,1),m)=sum(p);
                        end
                        for m=1:4
                            for n=1:16
                                p(n)=D(k,h(k,1),a(1,n))*D(k,h(k,3),a(2,n))*exp(-1/(2*sigma^2)*(abs(y(k)-S(k,h(k,2),m)-S(k,h(k,1),a(1,n))-S(k,h(k,3),a(2,n))))^2);  
                            end
                            Dtemp(k,h(k,2),m)=sum(p);
                        end
                        for m=1:4
                            for n=1:16
                                p(n)=D(k,h(k,1),a(1,n))*D(k,h(k,2),a(2,n))*exp(-1/(2*sigma^2)*(abs(y(k)-S(k,h(k,3),m)-S(k,h(k,1),a(1,n))-S(k,h(k,2),a(2,n))))^2);
                            end
                            Dtemp(k,h(k,3),m)=sum(p);
                        end
                    end
                    D=Dtemp;
                end        
                for j=1:J
                    De(1,j,:)=D(v(1,j),j,:).*D(v(2,j),j,:);
                    De(1,j,:)=De(1,j,:)/sum(De(1,j,:));
                    for m=1:4
                        if De(1,j,m)>0.9999
                            De(1,j,m)=0.9999;
                        elseif De(1,j,m)<0.0001
                            De(1,j,m)=0.0001;
                        end
                    end
                end
                for p2=1:J
                    Isoft_demod_output(2*p1-1:2*p1,p2)=[log((De(1,p2,2)+De(1,p2,1))/(De(1,p2,3)+De(1,p2,4))),log((De(1,p2,1)+De(1,p2,3))/(De(1,p2,2)+De(1,p2,4)))];
                end
%                 for j=1:J
%                     [~,receivedSignal(p1,j)]=max(De(1,j,:));
%                 end    
            end    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            toc
            %%%%%%%%%%% LDPC decoding %%%%%%%%%%%%%%%%
            for j=1:J
                soft_demod_output(:,j)=Isoft_demod_output(alpha(:,j),j);
                LdpcDecode(:,j) = step(hDec,soft_demod_output(:,j));
                ILdpcDecode(alpha(:,j),j)=LdpcDecode(:,j);
            end
            toc;
            for j=1:J
                for p1=1:LdpcLen/2
                    pp(p1,j,:)=[(1-1/(exp(ILdpcDecode(2*p1-1,j))+1))*(1-1/(exp(ILdpcDecode(2*p1,j))+1)),
                                (1-1/(exp(ILdpcDecode(2*p1-1,j))+1))*(1/(exp(ILdpcDecode(2*p1,j))+1)),
                                (1/(exp(ILdpcDecode(2*p1-1,j))+1))*(1-1/(exp(ILdpcDecode(2*p1,j))+1)),
                                (1/(exp(ILdpcDecode(2*p1-1,j))+1))*(1/(exp(ILdpcDecode(2*p1,j))+1))];
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        for p1=1:dataLen
            for j=1:J
                if LdpcDecode(p1,j)>0
                    receivedBits(p1,j)=0;
                else
                    receivedBits(p1,j)=1;
                end                  
            end
        end
        for j=1:J
            if data(:,j)==receivedBits(:,j)
                errBlc;
            else
                errBlc=errBlc+1;
            end
        end   
        if errBlc>=25
            flag=0;
        end
        errBit=errBit+sum(sum(abs(data-receivedBits)));        
        uerrBit(ii,:)=uerrBit(ii,:)+sum(abs(data-receivedBits),1);        
    end
    ber(ii)=errBit/J/dataLen/total
    bler(ii)=errBlc/J/total
    uber(ii,:)=uerrBit(ii,:)./dataLen/total
    uBber(ii)=min(uber(ii,:))
    toc;
end
toc;
save