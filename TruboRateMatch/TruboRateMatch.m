clear;
% clc;
R=1/2;
% 576,768,960,1152,1728,2304
rateMatchLen=576;
CRCdataLen=rateMatchLen*R;
frmLen = CRCdataLen-24;

rng default
EbNo= -1:0.5:10;
intrlvrIndices = randperm(CRCdataLen);

hTEnc = comm.TurboEncoder('TrellisStructure',poly2trellis(4, ...
    [13 15],13),'InterleaverIndices',intrlvrIndices);
hTDec = comm.TurboDecoder('TrellisStructure',poly2trellis(4, ...
    [13 15],13),'Algorithm','max','InterleaverIndices',intrlvrIndices, ...
    'NumIterations',8);
hMod = comm.BPSKModulator;
hGen = comm.CRCGenerator([24 23 18 17 14 11 10 7 6 5 4 3 1 0]);
hDetect = comm.CRCDetector([24 23 18 17 14 11 10 7 6 5 4 3 1 0]);
R=ceil((frmLen+4)/32);
C=32;
Ncb=CRCdataLen+4;
rv=0;
k0=R*(2*rv*(floor((Ncb+8*R-1)/(8*R)))+2);	
tic;
for ii=1:length(EbNo)
    SNR=EbNo(ii)
    noiseVar = 10^(SNR/10);
    hDemod = comm.BPSKDemodulator('DecisionMethod','Log-likelihood ratio', ...
    'Variance',noiseVar);
    % sigma=sqrt(1/2*noiseVar);
    hChan = comm.AWGNChannel('EbNo',SNR);
    hError = comm.ErrorRate;
    flag=1;
    frmIdx=0;
    errFrm=0;
    while flag
        frmIdx=frmIdx+1;
        data = randi([0 1],frmLen,1);
        CRCdata = step(hGen, data);
        encodedData = step(hTEnc,CRCdata);
        [systemBit,parityBit,iParityBit]=tiqu(encodedData);      
        rateMatchOut = rate_match( rateMatchLen,CRCdataLen,0,0,systemBit,parityBit,iParityBit,R,C,k0 );
        modSignal = step(hMod,rateMatchOut');
        receivedSignal = step(hChan,modSignal);
        demodSignal = step(hDemod,receivedSignal);
        [a b c]=DeRateMatch(demodSignal,CRCdataLen+4,k0,R,C,length(rateMatchOut),0,0);
        d=huanyuan(a,b,c);
        TurboDecode = step(hTDec,-1*d');
        receivedBits = step(hDetect, TurboDecode); 
        if receivedBits==data
            errFrm=errFrm;
        else
            errFrm=errFrm+1;
        end
        errorStats = step(hError,data,receivedBits);
        if errFrm>=50
            flag=0;
        end  
    end
    bler(ii)=errFrm/frmIdx
    ratio(ii)=errorStats(1)
    toc
end
toc;
semilogy(EbNo,ratio,'bo-');
grid on;
xlabel('EsNo');ylabel('BLER');
title('Length:Data=72 CRC=72+24 TurboRateMatch=192');