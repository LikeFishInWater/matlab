function Err = myTurboCodeSimu(K, Ncb, Rc, EbN0dB, Modulation, ChlType, Nits)
    

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
    intrlvrIndices = randperm(K);
    hTEnc = comm.TurboEncoder('TrellisStructure',poly2trellis(4, ...
    [13 15 17],13),'InterleaverIndices',intrlvrIndices);
    hTDec = comm.TurboDecoder('TrellisStructure',poly2trellis(4, ...
    [13 15 17],13),'InterleaverIndices',intrlvrIndices, ...
    'NumIterations',4);
    hDemod = comm.BPSKDemodulator('DecisionMethod','Log-likelihood ratio','Variance',Sigma2);
    for j = 1:Ncb
        j
        encodedData = step(hTEnc,InfoBit(:,j));              
        TxSymbol = 2*encodedData-1;
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
        SoftBits = step(hDemod,(RxSymbol.*conj(h))); 
        
        %Turbo decode
        DecodedBit = step(hTDec,-SoftBits);
        Err = Err + sum(xor(InfoBit(:,j),DecodedBit));
    end
end