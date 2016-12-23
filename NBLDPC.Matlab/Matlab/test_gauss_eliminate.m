clear; clc;


%% Initial
NonBinaryFileName = '200.300.regular.gf.4.txt';
qAry = 4;

% global MULTIPLY_TABLE;
% global INVERSE_TABLE;
NonBinary_LDPC_GenerateTable(qAry);

[H, N, M, Link_Var_to_Chk, Link_Chk_to_Var, qAry] = NonBinary_LDPC_LoadMatrixFile(NonBinaryFileName);
tic
[Eliminated_H, GaussEliminate_R, GaussEliminate_C, ExchangedOriginCol, ExchangedDestinCol]  = NonBinary_LDPC_GaussEliminate(H, qAry);
toc
spy(H); figure; spy(Eliminated_H);

EncodeInput = randi(qAry, 1, N-M) - 1;
EncodeOutput = NonBinary_LDPC_Encoder_GaussEliminate(EncodeInput, Eliminated_H, GaussEliminate_R, GaussEliminate_C, ExchangedOriginCol, ExchangedDestinCol, qAry);
isCorrect = NonBinary_LDPC_Encoder_Check(H, Link_Var_to_Chk, Link_Chk_to_Var, qAry, EncodeInput, EncodeOutput);