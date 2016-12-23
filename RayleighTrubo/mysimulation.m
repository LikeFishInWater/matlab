% Demonstration of Eb/N0 Vs BER for BPSK modulation scheme
clear;
%---------Input Fields------------------------

K = 2044; % Length of Code block
N = K; %Number of input bits
N0= (K+4)*3;
ModulType = 'BPSK';
ChlType = 'RAYLEIGH';
Rc=1/3; %Rc = code rate for a coded system. Since no coding is used Rc=1
% EbN0dB = -3.5:0.1:-2.7;
EbN0dB=-2:0.1:-1.5;

intrlvrIndices = randperm(K);
hTEnc = comm.TurboEncoder('TrellisStructure',poly2trellis(4,[13 15],13),'InterleaverIndices',intrlvrIndices);
hTDec = comm.TurboDecoder('TrellisStructure',poly2trellis(4,[13 15],13),'InterleaverIndices',intrlvrIndices,'NumIterations',4);
hMod = comm.BPSKModulator;
Err = zeros(1,length(EbN0dB)); %Place holder for BER values for each Eb/N0
Herr = zeros(1,length(EbN0dB));
BER_TC_AWGN = zeros(1,length(EbN0dB));
HBER_TC_AWGN = zeros(1,length(EbN0dB));
tic;
for i=1:length(EbN0dB)
    Rm = 1;
    flag=1;
    total=0;
    errFrm=0;
    EbN0 = 10.^(EbN0dB(i)/10); %Converting Eb/N0 dB value to linear scale
    Sigma2 = 1./(2*Rm*EbN0);
    noiseSigma = sqrt(Sigma2); %Standard deviation for AWGN Noise
    while flag
        total=total+1;
        InfoBit=randi([0 1],K,1); %Generating a uniformly distributed random 1s and 0s
        encodedData = step(hTEnc,InfoBit);    
        TxSymbol = step(hMod,encodedData);
        noise = noiseSigma*(randn(length(TxSymbol),1)+1i*randn(length(TxSymbol),1));
        switch ChlType
            case 'AWGN'
                h = 1;
            case 'RAYLEIGH'
                % h = 0.707*(randn(1)+1i*randn(1));
                h = 1/sqrt(pi/2)*(randn(length(TxSymbol),1) + 1i * randn(length(TxSymbol),1));
            otherwise
                disp('Not Support Channel Type');
                exit(-1);
        end            
        RxSymbol = h.*TxSymbol + noise;      
        %SoftBits = ((abs(RxSymbol+h)).^2-(abs(RxSymbol-h)).^2)/2/Sigma2;
        SoftBits = log(abs(RxSymbol-h)./abs(RxSymbol+h))-((abs(RxSymbol-h)).^2-(abs(RxSymbol+h)).^2)/2/Sigma2;
        DecodedBit = step(hTDec,-SoftBits);  
        for ii=1:length(SoftBits)
            if SoftBits(ii)>0
                hardDecode(ii)=0;
            else
                hardDecode(ii)=1;
            end
        end
        if InfoBit==DecodedBit
            errFrm=errFrm;
        else
            errFrm=errFrm+1;
        end
        if errFrm==50
            flag=0;
        end
        Err(i) = Err(i) + sum(xor(InfoBit,DecodedBit));
        Herr(i) = Herr(i) + sum(xor(encodedData,hardDecode')); 
    end
    BER_TC_AWGN(i) = Err(i)./N/total
    HBER_TC_AWGN(i) = Herr(i)./N0/total
    toc;
end
toc;

semilogy(EbN0dB,BER_TC_AWGN,'bs:','LineWidth',2);
hold on;
semilogy(EbN0dB,HBER_TC_AWGN,'rs:','LineWidth',2);
grid on;
legend('TC','UC');
