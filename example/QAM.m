clear;
M=4;
bitsPerSymbol=log(M)/log(2);

hMod = comm.RectangularQAMModulator('ModulationOrder',M,'BitInput',true,'NormalizationMethod','Average power');
hDemod = comm.RectangularQAMDemodulator('ModulationOrder',M,'BitOutput',true);

EbNo=1:10;
EsNo=EbNo+10*log10(bitsPerSymbol);
for ii=1:length(EbNo)
    snr=EsNo(ii);
    hAWGN = comm.AWGNChannel('NoiseMethod', ...
          'Signal to noise ratio (SNR)','SNR',snr);
    hError = comm.ErrorRate;
    flag=1;
    while flag
        dataIn = randi([0 1],60,1);
        txSig = step(hMod,dataIn);
        rxSig = step(hAWGN,txSig);
        dataOut = step(hDemod,rxSig);
        errorStats = step(hError,dataIn,dataOut);
        if errorStats(2)>1000
            flag=0;
        end
    end
    ber(ii)=errorStats(1)
end
semilogy(EbNo,ber);
grid on;

