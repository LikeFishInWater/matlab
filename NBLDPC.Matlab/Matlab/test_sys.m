clear; clc;


%% Initial
NonBinaryFileName = 'mn_504GF4.txt';
qAry = 4;
%Pre-caculate the gf(q) multiply and inverse table
NonBinary_LDPC_GenerateTable(qAry);

[H, N, M, Link_Var_to_Chk, Link_Chk_to_Var, qAry] = NonBinary_LDPC_LoadMatrixFile(NonBinaryFileName);
[Eliminated_H, GaussEliminate_R, GaussEliminate_C, ExchangedOriginCol, ExchangedDestinCol]  = NonBinary_LDPC_GaussEliminate(H, qAry);

%% Encode Random Msg
EncodeInput = randi(qAry, 1, N-M) - 1;
EncodeOutput = NonBinary_LDPC_Encoder_GaussEliminate(EncodeInput, Eliminated_H, GaussEliminate_R, GaussEliminate_C, ExchangedOriginCol, ExchangedDestinCol, qAry);
isCorrect = NonBinary_LDPC_Encoder_Check(H, Link_Var_to_Chk, Link_Chk_to_Var, qAry, EncodeInput, EncodeOutput);

%% Modulate
symbol_tx = zeros(1,N);
x0 = 1+1j;
x1 = 1-1j;
x2 = -1+1j;
x3 = -1-1j;
for i=1:N
    switch EncodeOutput(1,i)
        case 0,
            symbol_tx(1,i) = x0;
        case 1,
            symbol_tx(1,i) = x1;
        case 2,
            symbol_tx(1,i) = x2;
        case 3,
            symbol_tx(1,i) = x3;
    end
end

%% Channel
sigma = 0.5;
symbol_rx = symbol_tx + (randn(1,N) + 1j * randn(1,N)) * sigma;

%% Demodulate
L_ch0 = exp(-abs(symbol_rx - x0).^2/(sigma^2));
L_ch1 = exp(-abs(symbol_rx - x1).^2/(sigma^2));
L_ch2 = exp(-abs(symbol_rx - x2).^2/(sigma^2));
L_ch3 = exp(-abs(symbol_rx - x3).^2/(sigma^2));
L_ch = [log(L_ch1./L_ch0); log(L_ch2./L_ch0); log(L_ch3./L_ch0)];

%% Deocde with BP algorithm
MaxIter = 5;
tic;
DecodeOutput = NonBinary_LDPC_Decoder_BP(H, Link_Var_to_Chk, Link_Chk_to_Var, qAry, L_ch, MaxIter);
toc;

%% Caculate the error code symbols
errCount = 0;
for i=1:N
    if DecodeOutput(i) ~= EncodeOutput(i)
        errCount = errCount + 1;
    end
end
SER = errCount / N