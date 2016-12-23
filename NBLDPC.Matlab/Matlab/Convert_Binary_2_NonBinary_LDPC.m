clear; clc;


%% File Name
BinaryFileName = './matrix/mat_N380_J2K20_FSEMPEG_ShaBN20Z19R910.txt';
NonBinaryFileName = './matrix/mat_N380_J2K20_FSEMPEG_ShaBN20Z19R910.GF.64.txt';
qAry = 64;

%% Generate the GF(q) table
% NonBinary_LDPC_GenerateTable(qAry);
load GF.64.mat
global MULTIPLY_TABLE

%% Convert the binary to nonbinary
NonBinary_LDPC_ConvertMatrixFile(BinaryFileName, NonBinaryFileName, qAry);
% NonBinary_LDPC_ConvertMatrixFile_Optimized(BinaryFileName, NonBinaryFileName, qAry);

%% Load the converted nonbinary LDPC file, and plot the non-0 element
[H, N, M, Link_Var_to_Chk, Link_Chk_to_Var, qAry] = NonBinary_LDPC_LoadMatrixFile(NonBinaryFileName);
spy(H);