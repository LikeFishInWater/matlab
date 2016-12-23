clear;

dataLen=2000;
LdpcLen=dataLen*2;
BinaryFile='MN_4000.txt';
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
hDec = comm.LDPCDecoder('ParityCheckMatrix' , H,'OutputValue','Whole codeword','DecisionMethod','Hard decision','IterationTerminationCondition','Parity check satisfied','FinalParityChecksOutputPort',true,'MaximumIterationCount',50);
% hEnc = comm.LDPCEncoder;
% hDec = comm.LDPCDecoder;
hMod = comm.BPSKModulator;
hDemod = comm.BPSKDemodulator('DecisionMethod','Log-likelihood ratio');
EsNo=-3:0.2:-1;
tic;
for ii=1:length(EsNo)
    snr=EsNo(ii);
    sigma=sqrt(1/2*1/(10^(snr/10)));
    total=0;
    errBit=0;
    errFrm=0;
    flag=1;
    while(flag)
        total=total+1;
        x=randi([0 1],dataLen,1);
        y=step(hEnc,x);
        modSignal=step(hMod,y);
        yr=awgn(modSignal,snr);
        deModSignal=step(hDemod,yr);
        decode=step(hDec,deModSignal/(2*sigma^2));
        correct=1;
        for jj=1:dataLen
            if(decode(jj)==x(jj))
                errBit=errBit;
            else 
                errBit=errBit+1;
                correct=0;
            end
        end
        if correct==0
            errFrm=errFrm+1;
        end
        if (errFrm>30)
            flag=0;
        end
    end   
    ber(ii)=errBit/total/dataLen
    fer(ii)=errFrm/total
    toc;
end
semilogy(EsNo,ber);
xlabel('EsNo');
ylabel('ber');
grid on;