function Err = UncodedSystemSimu(N, EbN0dB, Modulation, ChlType)
    

    Rc = 1;
    
    switch Modulation
        case 'BPSK'
            Rm = 1;
        case 'QPSK'
            Rm = 2;
        case '16QAM'
            Rm = 4;
        case '64QAM'
            Rm = 6;
        case '256QAM'
            Rm = 8;
        otherwise
            disp('Not Support Modulation Type');
            exit(-1);
    end
            
    InfoBit=randi([0 1],N,1); %Generating a uniformly distributed random 1s and 0s
    EbN0 = 10.^(EbN0dB/10); %Converting Eb/N0 dB value to linear scale
    Sigma2 = 1./(2*Rm*Rc*EbN0);
    noiseSigma = sqrt(Sigma2); %Standard deviation for AWGN Noise
    Err = 0;
        
    %Modulated
    TxSymbol = lteSymbolModulate(InfoBit,Modulation);
    %Adding noise 
    noise = noiseSigma*(randn(length(TxSymbol),1)+1i*randn(length(TxSymbol),1));
    switch ChlType
        case 'AWGN'
            h = 1;
        case 'RAYLEIGH'
            h = 0.707*(randn(length(TxSymbol),1)+1i*randn(length(TxSymbol),1));
        otherwise
            disp('Not Support Channel Type');
            exit(-1);
    end            
    RxSymbol = h.*TxSymbol + noise;
    RxSymbol = RxSymbol.*conj(h);
        
    %demodulated
    Bits = lteSymbolDemodulate(RxSymbol,Modulation,'Hard');
    Err = Err + sum(xor(InfoBit,Bits));
    
end