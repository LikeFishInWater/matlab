clear;
clc;

%% Initial Simulation Parameters
MaxIter = 500;
Eb_N0_dB=[1.2:0.2:3];
CodeFileName = './matrix/L20R2A.K.512.GF.32.txt';

nsample = 1000000;
[qAry_cpp, TABLE_MULTIPLY_cpp, TABLE_ADD_cpp, TABLE_INVERSE_cpp, maxVarDegree_cpp, maxChkDegree_cpp, VarDegree_cpp, ChkDegree_cpp, ...
    H_cpp, Row_Link_Col_cpp, Col_Link_Row_cpp,Q_cpp, R_cpp, L_Post_cpp, LDR_Vector_cpp, L_SIGMA_cpp, L_RHO_cpp,CodeLen_cpp, ChkLen_cpp, MsgLen_cpp, ...
    Rate_cpp] = Initial_CPP(CodeFileName, './TABLE/Table.GF.32.txt');


%% Initial NBLDPC
NonBinaryFileName = CodeFileName;
qAry = 32;
load ./TABLE/GF.32.TABLE.mat
% NonBinary_LDPC_GenerateTable(qAry);
[H, N, M, Link_Var_to_Chk, Link_Chk_to_Var, qAry] = NonBinary_LDPC_LoadMatrixFile(NonBinaryFileName);
[Eliminated_H, GaussEliminate_R, GaussEliminate_C, ExchangedOriginCol, ExchangedDestinCol]  = NonBinary_LDPC_GaussEliminate(H, qAry);
P = [];
T = [];
P_sum = 0;
T_sum = 0;
for col=1:N
    if size(Link_Chk_to_Var{col},1) == 4
        P = [P,col];
        P_sum = P_sum + 1;
    else
        T = [T,col];
        T_sum = T_sum + 1;
    end
end
R = (N - M) / (N - P_sum)


%% Initial Modulator & Demodulator
ModType =1; %4qam
hMod= comm.RectangularQAMModulator('ModulationOrder',2^ModType,'BitInput',true,'NormalizationMethod',...
    'Average power');
hDemod = comm.RectangularQAMDemodulator('ModulationOrder',2^ModType,'BitOutput',true,'NormalizationMethod',...
    'Average power','DecisionMethod','Log-likelihood ratio');

%% Initial Others
%sum_min=0;
BER=zeros(1,length(Eb_N0_dB)) ;
FER=zeros(1,length(Eb_N0_dB)) ;
err_bit=zeros(nsample,length(Eb_N0_dB));
err_bit_sum=zeros(1,length(Eb_N0_dB)) ;
err_Frames=zeros(1,length(Eb_N0_dB)) ;
countCycles=zeros(1,length(Eb_N0_dB)) ;

Pr_simple=zeros(qAry,1);
Pr_LLR=zeros(qAry-1,N);

CodeBits = log2(qAry);
Modbits = ModType;

%% Simulation
for Eb_N0_i=1:length(Eb_N0_dB)
    sigma_n = 1 / sqrt(2 * R * ModType * 10^(Eb_N0_dB(Eb_N0_i) / 10));
    for ns=1:nsample 
        % NBLDPC Encoder
%         EncodeInput = randi(qAry, 1, N-M) - 1;
        EncodeInput = zeros(1, N-M);
        EncodeOutput = NonBinary_LDPC_Encoder_GaussEliminate(EncodeInput, Eliminated_H, GaussEliminate_R, GaussEliminate_C, ExchangedOriginCol, ExchangedDestinCol, qAry);
        % Modulation
        EncodeOutput_Bit = dec2bin(EncodeOutput, log2(qAry));
        [a, b] = size(EncodeOutput_Bit);
        EncodeOutput_Bit = str2num(reshape(EncodeOutput_Bit',a * b,1));
        ModSymbol = 1 - 2 * EncodeOutput_Bit;
        % Channel
        noise = randn(N * CodeBits  ,1) * sigma_n;
        Y = ModSymbol+ noise;
        % Demodulation
        Bit_LLR = 2 * (real(Y)) / (sigma_n^2);
        Bit_Pr_0 = exp(Bit_LLR) ./ (1 + exp(Bit_LLR));
        Bit_Pr_1 = 1 ./ (1 + exp(Bit_LLR));
        Bit_Pr_0 = reshape(Bit_Pr_0', CodeBits, N);
        Bit_Pr_1 = reshape(Bit_Pr_1', CodeBits, N);
        Code_Pr(1, :) = prod(Bit_Pr_0, 1);
        for x = 1:qAry-1
            code = flipdim(bitget(x, 1:CodeBits)' * ones(1, N), 1);
            Code_Pr(x + 1, :) = prod((1 - code) .* Bit_Pr_0 + code .* Bit_Pr_1, 1);
            Pr_LLR(x,:) = log(Code_Pr(x + 1, :) ./ Code_Pr(1, :));
        end
        Pr_LLR(:,P) = zeros(qAry-1, P_sum);
        [DecodeOutput_uint, DecodeOutput]= Decode_CPP(Pr_LLR, CodeLen_cpp, ChkLen_cpp, VarDegree_cpp, ChkDegree_cpp, ...
            Col_Link_Row_cpp, Row_Link_Col_cpp, qAry_cpp, TABLE_MULTIPLY_cpp, TABLE_ADD_cpp, TABLE_INVERSE_cpp, ...
            H_cpp, Q_cpp, R_cpp, L_Post_cpp, LDR_Vector_cpp, L_SIGMA_cpp, L_RHO_cpp, MaxIter);
%         biterr(DecodeOutput1, DecodeOutput)
        % Error Counting
       [err_bit(ns, Eb_N0_i), ~] = biterr(EncodeOutput(1,T), DecodeOutput(1,T), CodeBits);
        if err_bit(ns, Eb_N0_i)  ~= 0
            err_Frames(Eb_N0_i) = err_Frames(Eb_N0_i) + 1;
            err_bit_sum(Eb_N0_i) = err_bit(ns,Eb_N0_i) + err_bit_sum(Eb_N0_i) ;
        end
        countCycles(Eb_N0_i) = countCycles(Eb_N0_i) + 1;
        if mod(countCycles(Eb_N0_i),100) == 0
            disp(['GF:' num2str(qAry)   'Eb/N0: ' num2str(Eb_N0_dB(Eb_N0_i)) ' BER:' num2str(err_bit_sum(Eb_N0_i) / (countCycles(Eb_N0_i) * T_sum * CodeBits)) ...
                ' FER:' num2str(err_Frames(Eb_N0_i) / countCycles(Eb_N0_i))  ' errF:' num2str(err_Frames(Eb_N0_i))  '   Counts:'  num2str(countCycles(Eb_N0_i))]);
        end
        if err_Frames(Eb_N0_i) >=100 && countCycles(Eb_N0_i)  >= 100
            break;
        end
    end
    BER(Eb_N0_i) = err_bit_sum(Eb_N0_i) / (countCycles(Eb_N0_i) * T_sum * CodeBits);
    FER(Eb_N0_i) = err_Frames(Eb_N0_i) / countCycles(Eb_N0_i);
    disp(['*********************************' 'GF:' num2str(qAry)  ' snr: ' num2str(Eb_N0_dB(Eb_N0_i)) ' BER:' num2str(BER(Eb_N0_i)) ...
                ' FER:' num2str(FER(Eb_N0_i))    '   Counts:'  num2str(countCycles(Eb_N0_i))]);
    save(mfilename)
end
save(mfilename,'Eb_N0_dB','BER','FER','countCycles','err_Frames','err_bit_sum');