ModM = 2;
SNRindB = 3.8;
CH_mode = 1; %0:AWGN;1:Rayleigh
r = 1/2;
% load('H_short_1_3.mat');
BinaryFile = 'MN_4000.txt';
fid=fopen(BinaryFile);
line=str2num(fgetl(fid));
Nn=line(1);
Mm=line(2); 
Rr = (Nn - Mm) / Nn;
fgetl(fid);fgetl(fid);fgetl(fid);
H_row_id = [];
H_col_id = [];
for col = 1:Nn
    row_id = str2num(fgetl(fid));
    row_id(row_id == 0) = [];
    H_row_id = [H_row_id, row_id];
    col_id = col * ones(1,length(row_id));
    H_col_id=[H_col_id, col_id];
end
fclose(fid);
H = sparse(H_row_id, H_col_id, ones(1,length(H_row_id)));
H = (H >= 1);
clear H_row_id;clear H_col_id;

filename = strcat('link_short_13_','.xls');
filename_temp = strcat('link_short_13_','temp','.xls');

N = 4000;
K = N*r;

MaxFEN = 25;  %������MaxFEN��������
FECNum = 1000;

modObj = modem.qammod('M',2^ModM,'SymbolOrder', 'Gray','InputType','Bit');
codeword = de2bi([0:2^ModM-1]);
CON = modulate(modObj, reshape(codeword.',1,[]).').';

DeNum = 1;

for idx = 1:length(SNRindB)
    SNRdB = SNRindB(idx);
    result_M = zeros(1,2);
    BlockNum = 0;
    FEN1 = 0;
    FEN2 = 0;
    result_temp = [];
    
    while FEN1 < MaxFEN
        result_inpar = zeros(1,2);
        [FEN1_temp BEN1_temp] = LDPC_NoBCH_DeNum_NUC_f(r,N,ModM,SNRdB,H,CON,CH_mode,DeNum);
        result_inpar = [FEN1_temp BEN1_temp];
        BlockNum = BlockNum + DeNum;
        result_M = result_M + result_inpar;
        FEN1 = result_M(:,1);
        BEN1 = result_M(:,2);
        result_temp = [BlockNum result_M];
        disp(result_temp);
    end
    result(idx,:)= [SNRdB BlockNum FEN1 BEN1 FEN1/BlockNum  BEN1/BlockNum/(N*r)]
    BEN1/BlockNum/(N*r)
end

