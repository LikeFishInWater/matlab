clear;
EsNo= -4:0.2:-3;
frmLen = 256;
rng default
noiseVar = 10.^(-EsNo./10);
intrlvrIndices = randperm(frmLen);
hTEnc = comm.TurboEncoder('TrellisStructure',poly2trellis(4, [13 15],13),'InterleaverIndices',intrlvrIndices);%%rate=1/3;
hTDec = comm.TurboDecoder('TrellisStructure',poly2trellis(4, [13 15],13),'InterleaverIndices',intrlvrIndices, 'NumIterations',4);
hMod = comm.BPSKModulator;
total=zeros(1,length(EsNo));
errFrm=zeros(1,length(EsNo));
errBit=zeros(1,length(EsNo));
for i=1:length(EsNo)
    flag=1;
    hChan = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (Es/No)','EsNo',EsNo(i));
    hDemod = comm.BPSKDemodulator('DecisionMethod','Log-likelihood ratio','Variance',noiseVar(i));
    while flag
        total(i)=total(i)+1;
        data = randi([0 1],frmLen,1);
        encodedData = step(hTEnc,data);
        modSignal = step(hMod,encodedData);
        receivedSignal = step(hChan,modSignal);
        demodSignal = step(hDemod,receivedSignal);
        receivedBits = step(hTDec,-demodSignal);
        if receivedBits == data
            errFrm(i)=errFrm(i);
        else 
            errFrm(i)=errFrm(i)+1;
        end
        if errFrm(i)==50
            flag=0;
        end
        errBit(i)=errBit(i)+nnz(receivedBits-data);
    end
    ratio_fe(i)=errFrm(i)/total(i)
    ratio_be(i)=errBit(i)/total(i)/frmLen
end

semilogy(EsNo,ratio_be,'r*-');
grid on;
xlabel('EsNo');
ylabel('BER');


