close all;
clear all;

randn('state',sum(100*clock));

%fileLoad='wu_matrix.mat';

fileLoad='davey1_data_3_5_504';%(j=3,k=5,CodeRate=1/3)

%fileLoad='davey1_data_3_9_2016.mat';
%fileLoad='davey1_data2016_1_2_63.mat';
%fileLoad='davey1_data_1_2_63.mat';%'gallager1_data.mat';% 'davey1_data.mat' ;
loadData=load(fileLoad);

rankH=loadData.rankA;
n=loadData.N0;
k=loadData.K0;
m=loadData.M0;
CodeRate=loadData.CodeRate;

if rankH~=m
     disp('rank of H should equal to the number rows of H.');    
     return; 
end    

H=loadData.H;
GT=loadData.GT;
clear loadData; %Free memory used by loadData



[childrenList,parentsList,children,parents]=RelationGraph(H);%Get the relation graph of the given LDPC

%Allocate map memory for the matrix R and Q (each element in the matrix has two status: 0,1)
[rows,columns]=size(parentsList);
R_map=zeros(rows,columns,2); %Allocate memory for the matrix R 
rowsRmap=rows;
colRmap=columns;


[rows,columns]=size(childrenList);
Q_map=zeros(rows,columns,2); %Allocate memory for the matrix Q
rowsQmap=rows;
colQmap=columns;

F0=zeros(n,2);%F0(:,1) for a=0;F0(:,2) for a=1

tentativeDecoder=zeros(n,2);

EbN0=4; %EbN0 dB

TestLength=10000000;

maxIterativeNumber=100;

A=1.0;

Eb=A*A/CodeRate;

N0=Eb*10^(-EbN0/10);
NoiseVar=N0/2;
NoiseRoot=sqrt(NoiseVar);

%Allocate memory for the variables
x=zeros(n,1); % Bipolar signal from set (+1,-1) to be transmitted
r=zeros(n,1); %Received signal polluted by noise
NoiseVector=zeros(n,1); %Noise vector
s_in=zeros(k,1); %Source binary sequence (sequence 0,1)
t_in=zeros(n,1); %Coded sequence ( sequence 0,1)

LLR=zeros(n,1);
valueLLR=zeros(n,1);
yDe=zeros(n,1);
yDe_maxSNR=zeros(n,1);

NoiseVector=randn(n,1)*NoiseRoot;
s_in=(randn(k,1)>0); %source information sequence with length K0
t_in=mod(GT*s_in,2);
x=(1-2*t_in);
r=A*x+NoiseVector;

%Estimate A and Noise Variance
r1=abs(r);
A1=mean(r1);
NoiseVar1=var(r1);

BlockCount=0;
IgnoreCount=1;

SourceCount=0;
SourceBlockLength=k;
ErrorCount=0;

ErrorCount1=0;

historyErrorFlag=1;
while(1) % start coding and decoding, stop when TestLength is achieved
    
    %******* Coding ********
    s_in=(randn(k,1)>0); %source information sequence with length k
    
    t_in=mod(GT*s_in,2);
    x=(1-2*t_in);
    %******** AWGN Channel ********
    NoiseVector=randn(n,1)*NoiseRoot;
    r=A*x+NoiseVector;
    
    %**************** Decoding ********
    A1=mean(abs(r));
    NoiseVar1=var(abs(r));
    
    valueFac=2*A1/NoiseVar1; %Initialize LLR from channel
    
    %valueFac=2*A/NoiseVar;
    
    LLR=valueFac*r;
    
    valueLLR=exp(LLR);
    F0(:,1)=valueLLR./(1+valueLLR);
    F0(:,2)=1-F0(:,1); %1./(1+valueLLR);
    
    
    %Initialize the matrix Q(i,j)  (Q_map)
    for u=1:1:rowsQmap
        for v=1:1:colQmap
            if(childrenList(u,v)>0)
                for a=1:1:2
                    Q_map(u,v,a)=F0(v,a);  
                end
            end
        end
    end
        

    Pr_x=zeros(max(parents),2); %the buffer for extrinsic information Pr(xj=0) and Pr(xj=1) (j=1,2,...)
    R_x=zeros(max(children),2);
    R_x_inv=zeros(max(children),2);
    F_j=zeros(max(parents),2);
    B_j=zeros(max(parents),2);

    maxSNR_LLR=0;
    errorFlag=1;
    %******** Iterative unit (start)****************
    for iterativeCount=1:1:maxIterativeNumber
        
        %Update the matrix R
        for u=1:1:rowsRmap %row u-th presents the u-th check equation
            L=parents(u);
            if(L==0) %current row doen't contain any child or parent (this will not occur for m equals to rank of H)
                continue;
            end
            for v=1:1:L
                i=u;%transfer to the locations of H
                j=parentsList(u,v);
                rowQ=sum(H(1:i,j));%transfer to the locations of Q_map
                colQ=j;
                for a=1:1:2
                    Pr_x(v,a)=Q_map(rowQ,colQ,a);
                end
            end
            %Calculate forward parameters
            F_j(1,1)=1;%Intialize F_j 
            F_j(1,2)=0;
            for j=2:1:L
                for s=1:1:2
                    F_j(j,s)=0;
                    for s1=1:1:2
                        tempStatus=(s-1)+(s1-1);
                        tempStatus=mod(tempStatus,2);
                        tempStatus=tempStatus+1;
                        F_j(j,s)=F_j(j,s)+Pr_x(j-1,tempStatus)*F_j(j-1,s1);  
                    end
                end
            end
            %Calculate backward parameters
            B_j(L,1)=1;%Intialize B_j
            B_j(L,2)=0;
            for j=(L-1):(-1):1
                for t=1:1:2
                    B_j(j,t)=0;
                    for t1=1:1:2
                         tempStatus=(t-1)+(t1-1);
                         tempStatus=mod(tempStatus,2);
                         tempStatus=tempStatus+1;
                         B_j(j,t)=B_j(j,t)+Pr_x(j+1,tempStatus)*B_j(j+1,t1);
                    end
                end
            end
            
            for j=1:1:L %update the matrix R
                for a=1:1:2
                    tempValue=0;
                    for s=1:1:2
                        tempStatus=mod((a-1)+(s-1),2);
                        tempStatus=tempStatus+1;
                        tempValue=tempValue+F_j(j,s)*B_j(j,tempStatus);    
                    end
                    R_map(u,j,a)=tempValue;
                end
            end
            
                   
        end

        
        %Update the matrix Q
        for v=1:1:n % (w,v) is the location in the Q_map
            for w=1:1:children(v)
                j=v; %transfer to the location of H  (i,j)
                i=childrenList(w,v);    
                rowR=i;
                colR=sum(H(i,1:j));
                for a=1:1:2
                    R_x(w,a)=R_map(rowR,colR,a);
                end
            end
            R_x_inv=1./(R_x+eps);
            productR=[1,1];
            for a=1:1:2
                for w=1:1:children(v)
                    productR(a)=productR(a)*R_x(w,a);    
                end
            end
            
            for w=1:1:children(v)%update the matrix Q
                tempValue=0;
                for a=1:1:2
                    Q_map(w,v,a)=F0(v,a)*productR(a)*R_x_inv(w,a);    
                    tempValue=tempValue+Q_map(w,v,a);    
                end
                tempValue=tempValue+eps;%divide 0 error process
                for a=1:1:2
                    Q_map(w,v,a)=Q_map(w,v,a)/tempValue;    
                end
            end
        end
                
        %Tentative decoding
        for v=1:1:n
            for a=1:1:2
                tentativeDecoder(v,a)=F0(v,a);
            end
            for w=1:1:children(v)    
                j=v;
                i=childrenList(w,v);
                rowR=i;
                colR=sum(H(i,1:j));
                for a=1:1:2
                    tentativeDecoder(v,a)=tentativeDecoder(v,a)*R_map(rowR,colR,a);
                end
            end
        end
        LLR=tentativeDecoder(:,1)./(tentativeDecoder(:,2)+eps);
        LLR=log(LLR);
        
        meanLLR=mean(abs(LLR));
        varLLR=var(abs(LLR));
        snrLLR=(meanLLR*meanLLR)/varLLR
        iterativeCount
        
        %Tentative Decoding
        yDe=(tentativeDecoder(:,1)<=tentativeDecoder(:,2));
        zDe=mod(H*yDe,2);
        valueDe=sum(zDe)
        if(valueDe==0)
            disp(iterativeCount)
            errorFlag=0;
            historyErrorFlag=0;
            break;
        end
        if(isnan(snrLLR) & (iterativeCount==1)) 
            disp('snrLLR is NaN');
            break;
        end
        if(snrLLR>maxSNR_LLR)
            maxSNR_LLR=snrLLR;
            yDe_maxSNR=yDe;
        end
        if(iterativeCount==maxIterativeNumber)
            yDe=yDe_maxSNR;    
        end
        
        
        
        
    end
    %***************Iterative unit (end) ***********
    
    
    
    
    
    %**** Re-estimate A1 and noiseVar1 ****
%     if(0==errorFlag)
%         r1=r.*(1-2*yDe);
%         A1=mean(r1)            
%         NoiseVar1=var(r1)
%     end
%     if(1==historyErrorFlag)
%         r1=r.*(1-2*yDe);
%         A1=mean(r1)            
%         NoiseVar1=var(r1)   
%     end
    
    %**** end of re-estimate A1 and NoiseVar1
    
    decoderOutput=yDe(1:k);
    decoderDirect=(r(1:k)<=0);

    BlockCount=BlockCount+1;
    iterativeTimes(BlockCount)=iterativeCount;
    if(BlockCount>IgnoreCount)
        errorVector=xor(s_in,decoderOutput);
        errornum=sum(errorVector);
        ErrorCount=ErrorCount+errornum;
        %if ((errornum>0) & (zDe==0) )
        %    disp('Low weight code encounted...');
        %    break;
        %end
        
        
        
        errorVector1=xor(s_in,decoderDirect);
        errornum1=sum(errorVector1);
        ErrorCount1=ErrorCount1+errornum1;
        
        
        
        SourceCount=SourceCount+SourceBlockLength; 
        
        BER_DirectOutput=ErrorCount1/SourceCount;
        BER=ErrorCount/SourceCount;
                
        BER_DirectOutput
        ErrorCount1
        
        BER
        ErrorCount
        SourceCount
        CurrentBlockIndex=BlockCount-IgnoreCount
    end
    
    if(SourceCount>=TestLength)
        break;
    end
end


