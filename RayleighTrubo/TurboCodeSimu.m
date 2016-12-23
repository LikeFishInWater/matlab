function Err = TurboCodeSimu(K, Ncb, Rc, EbN0dB, Modulation, ChlType, Nits)
    

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
            
    InfoBit=randi([0 1],K,Ncb); %Generating a uniformly distributed random 1s and 0s
    EbN0 = 10.^(EbN0dB/10); %Converting Eb/N0 dB value to linear scale
    Sigma2 = 1./(2*Rm*Rc*EbN0);
    noiseSigma = sqrt(Sigma2); %Standard deviation for AWGN Noise
    Err = 0;
    for j = 1:Ncb
        EncodedBit = lteTurboEncode(InfoBit(:,j));%Coded bits
                
        %bpskModulated = 2*EncodedBit-1; %Mapping 0-&gt;-1 and 1-&gt;1
        TxSymbol = lteSymbolModulate(EncodedBit,Modulation);
        %Adding noise 
        noise = noiseSigma*(randn(length(TxSymbol),1)+1i*randn(length(TxSymbol),1));
        switch ChlType
            case 'AWGN'
                h = 1;
            case 'RAYLEIGH'
                h = 0.707*(randn(1)+1i*randn(1));
            otherwise
                disp('Not Support Channel Type');
                exit(-1);
        end            
        RxSymbol = h.*TxSymbol + noise;
        RxSymbol = RxSymbol.*conj(h);
        ChlGain = conj(h)*h;
        
        %Soft demodulated
        SoftBits = 1.414/ChlGain/Sigma2*lteSymbolDemodulate(RxSymbol,Modulation,'Soft');
        
        %Turbo decode
        DecodedBit = lteTurboDecode(SoftBits,Nits);
        Err = Err + sum(xor(InfoBit(:,j),DecodedBit));
    end
end