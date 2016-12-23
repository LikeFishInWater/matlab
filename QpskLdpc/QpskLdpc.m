clear;
errBlcNum=50;
LdpcIte=30;
dataLen=972;
LdpcLen=1944;
EsNo=1:0.2:2;

BinaryFile='M_1944_1_2.txt';
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
hDec = comm.LDPCDecoder('ParityCheckMatrix' , H,'OutputValue','Whole codeword','IterationTerminationCondition','Parity check satisfied','FinalParityChecksOutputPort',true,'MaximumIterationCount',LdpcIte);


hMod = comm.QPSKModulator('PhaseOffset',pi/4,'BitInput',true);
%Create an error rate calculator
for ii=1:length(EsNo)
    snr=EsNo(ii)
    sigma=sqrt(1/2*10^(-snr/10));
    hAWGN = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)','SNR',snr);
    hDemod = comm.QPSKDemodulator('PhaseOffset',pi/4,'BitOutput',true,'DecisionMethod','Log-likelihood ratio','Variance',sigma);
    flag=1;
    total=0;
    errBit=0;
    errBlc=0;
    while flag
        total=total+1;
        data = randi([0 1],dataLen,1);
        LdpcData=step(hEnc,data);
        modSignal = step(hMod, LdpcData);
        % noisySignal = step(hAWGN, modSignal);
        noisySignal=awgn(modSignal,snr);
        SoftDemodOutput = step(hDemod, noisySignal);
        [receivedData,parity]=step(hDec,SoftDemodOutput);
        if receivedData==LdpcData
            errBlc=errBlc;
        else 
            errBlc=errBlc+1;
        end
        errBit=errBit+sum(sum(abs(receivedData-LdpcData)));        
        if errBlc==errBlcNum
            flag=0;
        end
    end
    ber(ii)=errBit/LdpcLen/total
    bler(ii)=errBlc/total
end 