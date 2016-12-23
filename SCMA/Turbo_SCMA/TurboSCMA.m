clear;
%%%%%%%%%% scma module %%%%%%%%%%%%%
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
load Sv1.mat;
h=zeros(4,3);
for m=1:4        
    h(m,:)=find(F(m,:));
end
v=zeros(2,6);
for m=1:6
	v(:,m)=find(F(:,m));
end
EsNo=2.5:0.2:7;

%%%%%%%%%%% turbo module 


intrlvrIndices = randperm(2044);
hTEnc = comm.TurboEncoder('TrellisStructure',poly2trellis(4,[13 15],13),'InterleaverIndices',intrlvrIndices);
hTDec = comm.TurboDecoder('TrellisStructure',poly2trellis(4,[13 15],13),'InterleaverIndices',intrlvrIndices,'NumIterations',4);

ModBits = 2;
hMod = comm.RectangularQAMModulator('ModulationOrder',2^ModBits,'BitInput',true,'NormalizationMethod', 'Average power');


%%%%%%%%%%%%%%%%%
tic;
r=2/5;
dataLen=2044;
TurboLen=(dataLen+4)*3;
rateMatchLen=dataLen/r;
data=zeros(dataLen,J);
encodedData=zeros(TurboLen,J);
rateMatchOut=zeros(rateMatchLen,J);
IrateMatchOut=zeros(rateMatchLen,J);
modSignal=zeros(rateMatchLen/2,J);
uerrBit=zeros(length(EsNo),J);
TurboDecode=zeros(dataLen,J);

R=ceil((dataLen+4)/32);
C=32;
Ncb=dataLen;
rv=0;
k0=R*(2*rv*(floor((Ncb+8*R-1)/(8*R)))+2);	
itera=8;
for j=1:J
    [temp, alpha(:,j)] = sort(rand(1,rateMatchLen));  
end
a=zeros(2,16);
for i=1:16
    a(1,i)=mod(i-1,4);
    a(2,i)=(i-1-a(1,i))/4;
end
a=a+1;
Dinit=zeros(K,J,M);
for ii=1:length(EsNo)
    snr=EsNo(ii)    
    flag=1;
    total=0;
    errBlc=0;
    errBit=0;
    %while flag 
    while flag
        pp=1/M*(zeros(rateMatchLen/2,J,M)+1);  
        total=total+1;
        %%%%%%%%%% Turbo coding and modulation module %%%%%%%%%%%
        for j=1:J
            data(:,j) = randi([0 1], dataLen, 1);
            encodedData(:,j) = step(hTEnc,data(:,j));
            [systemBit,parityBit,iParityBit]=tiqu(encodedData(:,j));  
            rateMatchOut(:,j) = rate_match( rateMatchLen,dataLen,0,0,systemBit,parityBit,iParityBit,R,C,k0 );
            IrateMatchOut(alpha(:,j),j)=rateMatchOut(:,j);
            for p1=1:rateMatchLen/2
                modSignal(p1,j)=bin2dec(sprintf('%d',IrateMatchOut(p1*2-1,j),IrateMatchOut(p1*2,j)))+1;
            end
        end       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% transmition and receiving %%%%%%%%%%%
        for p1=1:rateMatchLen/2
            x=zeros(K,1);
            for q1=1:J
                x=x+S(:,q1,modSignal(p1,q1));
            end
            yt(:,p1)=awgn(x,snr);
            sigma=sqrt(1/2*10^(-snr/10));
        end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% SCMA coding and decoding module %%%%%%%%%%%%
        for p1=1:rateMatchLen/2
            y=yt(:,p1);
            sigma=sqrt(1/2*10^(-snr/10));
            D=zeros(K,J,M);
            for m=1:4
                for n=1:3
                    D(m,h(m,n),:)=1/4*[1 1 1 1];
                end
            end
            T=8;
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

            De=zeros(1,J,M);       
            for j=1:J
                De(1,j,:)=D(v(1,j),j,:).*D(v(2,j),j,:).*pp(p1,j,:);
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
            for j=1:J
                [~,receivedSignal(p1,j)]=max(De(1,j,:));
            end     
        end    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%% Turbo decoding %%%%%%%%%%%%%%%%

        for j=1:J
            soft_demod_output(:,j)=Isoft_demod_output(alpha(:,j),j);
            [a0 b c]=DeRateMatch(soft_demod_output(:,j),dataLen+4,k0,R,C,length(rateMatchOut),0,0);
            d=huanyuan(a0,b,c);
            TurboDecode(:,j) = step(hTDec,-d');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j=1:J
            if data(:,j)==TurboDecode(:,j)
                errBlc;
            else
                errBlc=errBlc+1;
            end
        end   
        if errBlc>=25
            flag=0;
        end
        errBit=errBit+sum(sum(abs(data-TurboDecode)));        
        uerrBit(ii,:)=uerrBit(ii,:)+sum(abs(data-TurboDecode),1);        
    end
    ber(ii)=errBit/J/dataLen/total
    bler(ii)=errBlc/J/total
    uber(ii,:)=uerrBit(ii,:)./dataLen/total
    uBber(ii)=min(uber(ii,:))
    toc;
end
toc;
