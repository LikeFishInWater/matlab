clear;
hMod = comm.BPSKModulator;
hDemod = comm.BPSKDemodulator;




EsNo=0:10;
for ii=1:length(EsNo)
    snr=EsNo(ii);
    hAWGN = comm.AWGNChannel('NoiseMethod', ...
    'Signal to noise ratio (SNR)','SNR',snr);
    hError = comm.ErrorRate;
    flag=1;
    while(flag)
        data = randi([0 1],500,1);                    % Generate data
        modSignal = step(hMod,data);                 % Modulate
        noisySignal = step(hAWGN,modSignal);         % Pass through AWGN
        receivedData = step(hDemod,noisySignal);     % Demodulate
        errorStats = step(hError,data,receivedData); % Collect error stats
        if errorStats(2)>100
            flag=0;
        end
    end
    ber(ii)=errorStats(1)
end
semilogy(EsNo,ber)
grid on;
save bpsk