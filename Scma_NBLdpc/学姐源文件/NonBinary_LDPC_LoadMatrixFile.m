%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NAME:        NonBinary_LDPC_LoadMatrixFile
% PURPOSE:     Read the Nonbinary LDPC code file
%
% Input:
% FileName: NonBinary LDPC code file name
%
% Output:
% H: GF(qAry) NonBinary LDPC check Matrix
% R: R is the indexes of variable nodes connted to Chk_j
% C: C is the indexes of check nodes connted to Var_i
% qAry: GF(qAry)
% N: Code Length
% M: Check Length
% 
% AUTHOR:       Xiaoshi
% DATE:         2014.12.04
% VERSION:      v1.0
% REVISED BY:   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Load Parity Check Matrix file, transfer it into sparse matrix H
function [H, N, M, R, C, qAry] = NonBinary_LDPC_LoadMatrixFile(FileName)
% fid=fopen('mn_504GF4.txt');
fid=fopen(FileName);
line=str2num(fgetl(fid));
N = line(1); M = line(2); qAry = line(3);
fgetl(fid);fgetl(fid);fgetl(fid);

H_row_id = [];
H_col_id = [];
H_value = [];
for col = 1:N
    line = str2num(fgetl(fid));
    line(line == 0) = [];
    for d=1:size(line,2)
        if mod(d,2) == 1
            H_row_id =[H_row_id, line(d)];
        else
            H_value = [H_value, line(d)];
        end
    end
    col_id = col * ones(1,length(line)/2);
    H_col_id=[H_col_id, col_id];
end
fclose(fid);
H = sparse(H_row_id, H_col_id, H_value);

R=cell(M,1);
C=cell(1,N);
for j=1:M
    i_ext=find(H_row_id(1,:)==j);
    R{j}=H_col_id(i_ext);
end
for i=1:N
    j_ext=find(H_col_id(1,:)==i);
    C{i}=H_row_id(j_ext)';
end
