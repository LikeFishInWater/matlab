clear;
M=16;
bitPerSym=log(M)/log(2);    
% EbNo=0:20;
% EsNo=EbNo+10*log10(bitPerSym);
 EsNo=11:20;
hMod = comm.PSKModulator(M, 'PhaseOffset',pi/16,'BitInput',true);
hDemod = comm.PSKDemodulator(M, 'PhaseOffset',pi/16,'BitOutput',true);
for ii=1:length(EsNo)
    snr=EsNo(ii);
    total=0;
    flag=1;
    hAWGN = comm.AWGNChannel('NoiseMethod', ...
          'Signal to noise ratio (SNR)','SNR',snr);
    hError = comm.ErrorRate;
    while flag
        data = randi([0 1],60,1);
        modSignal = step(hMod, data);
        noisySignal = step(hAWGN, modSignal);
        receivedData = step(hDemod, noisySignal);
        errorStats = step(hError, data, receivedData);
        if errorStats(2)>50
            flag=0;
        end
    end
    % errBit=errorStats(2)
    BER(ii)=errorStats(1)
end
semilogy(EbNo,BER);
xlabel('EbNo(dB)');ylabel('BER');
grid on;
    
    
   