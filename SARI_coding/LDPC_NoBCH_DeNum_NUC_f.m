%输入星座图
function [FEN1 BEN1] = LDPC_NoBCH_DeNum_NUC_f(r,N,ModM,SNRindB,H,CON,CH_mode,DeNum)
%确定block数的链路仿真，block数由DeNum数指定
%只跑一个SNR点
%为了加快仿真速度，没有BCH编码

K = N*r;
M = N -K;

rand('state',sum(100*clock));

% Construct a default LDPC encoder object
% enc = fec.ldpcenc(H);
enc = comm.LDPCEncoder(H);
% Construct a companion LDPC decoder object
dec = comm.LDPCDecoder('ParityCheckMatrix',H,'DecisionMethod','Hard decision','OutputValue','Information part','IterationTerminationCondition','Parity check satisfied');
% dec = fec.ldpcdec(H);
% 
% dec.DecisionType = 'Hard decision';
% dec.OutputFormat = 'Information part';
% dec.NumIterations = 50;
% % Stop if all parity-checks are satisfied
% dec.DoParityChecks = 'Yes';


EbN0 = SNRindB - 10*log10(r*ModM);

SNRdB = SNRindB;
sigma = sqrt(10^(-SNRdB/10));

FEN1 = 0;
BEN1 = 0;
FEN2 = 0;
BEN2 = 0;

for BlockNum = 1: DeNum

    % msg = randint(1,K,2);
    msg = randi([0 1],K,1) ;
    % codeword = encode(enc,msg);
    codeword = step(enc,msg);
    
    % Modulate the signal
    modulatedsig = CON(bi2de(fliplr(reshape(codeword.', ModM, N/ModM).'))+1).';    
    
    % Transmit signal through channel
    
    if CH_mode == 0
        channel = ones(1,length(modulatedsig));
        N0 = sigma^2*ones(1,length(modulatedsig)).'*mean(abs(modulatedsig).^2);
    else
        channel = sqrt(0.5) * randn(1,length(modulatedsig)) + sqrt(0.5) * j * randn(1,length(modulatedsig));
        N0 = sigma^2./abs(channel).^2.'*mean(abs(modulatedsig).^2);
    end
    
    receivedsig = awgn(modulatedsig.*channel.', SNRdB, 'measured').';%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%?????????????????????????????????????????????
      
    llr = QAM_demap_NUC(receivedsig.'./channel.',N0,CON);
    
    % llr = llr.';
    
    % decodedmsgldpc = decode(dec, llr);%LDPC码译码
    decodedmsgldpc = step(dec, llr);%LDPC码译码
        
    if nnz(decodedmsgldpc-msg)>0
        FEN1 = FEN1+1;
        BEN1 = BEN1 + nnz(decodedmsgldpc-msg);
    end
    
    %         if mod(BlockNum,10) ==0
    %             fprintf('BlockNum = %d,  FENldpc = %d,  BENldpc = %d     iter_num = %d\n',BlockNum,FEN1,BEN1,dec.ActualNumIterations);
    %             fprintf('BlockNum = %d,  FENbch = %d,  BENbch = %d     iter_num = %d\n',BlockNum,FEN2,BEN2,dec.ActualNumIterations);
    %         end
end

