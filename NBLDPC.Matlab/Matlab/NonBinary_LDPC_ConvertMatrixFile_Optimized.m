%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NAME:        NonBinary_LDPC_ConvertMatrixFile
% PURPOSE:     Convert the Binary LDPC check matrix into Nonbinary by
% random replacement of the non-0 element && conver it into BiDirection
%
% Input:
% BinaryFileName: Binary LDPC code file name
% NonBinaryFileName: NonBinary LDPC code file name
%
% AUTHOR:       Xiaoshi
% DATE:         2014.12.04
% VERSION:      v1.0
% REVISED BY:   Xiaoshi 2015.03.17 remove the 'not used'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function NonBinary_LDPC_ConvertMatrixFile_Optimized(BinaryFileName, NonBinaryFileName, qAry)

global MULTIPLY_TABLE

%% Load Binary LDPC File
fid=fopen(BinaryFileName);
line=str2num(fgetl(fid));
N=line(1); M = line(2);
maxDegree = str2num(fgetl(fid));
colWei = str2num(fgetl(fid));
rowWei = str2num(fgetl(fid));

H_row_id = [];
H_col_id = [];
for col = 1:N
    row_id = str2num(fgetl(fid));
    row_id(row_id == 0) = [];
    H_row_id = [H_row_id, row_id];
    col_id = col * ones(1,length(row_id));
    H_col_id=[H_col_id, col_id];
end
fclose(fid);

H = sparse(H_row_id, H_col_id, ones(1,length(H_row_id)));
% H = (H >= 1);
clear H_row_id;clear H_col_id;

%%Convert the Binary into NonBinary LDPC file
%conver it into BiDirection
% for i=1:M
%     if i == 1
%         H(i,N-M+1:N) = [1 zeros(1,M-1)];
%     else
%         H(i,N-M+1:N) = [zeros(1,i-2) ones(1,2) zeros(1,M-i)];
%     end
% end

%% replace the non-0 element with random variable 1~q-1
dc4_tuple_gf32 = [31, 17, 5 ,1;
    12, 31, 5, 1;
    30, 31, 5, 1;
    30, 31, 10, 1;
    24, 29, 20, 1;
    15, 31, 20, 1;
    21, 29, 20, 1
    ];

dc5_tuple_gf16 = [14 5 3 1 1; 14 11 3 1 1;
14 5 3 2 1; 14 6 5 2 1; 14 10 6 2 1; 15 10 6 2 1;
14 11 6 2 1; 14 5 3 3 1; 15 5 3 3 1; 15 5 4 3 1;
15 5 5 3 1; 15 6 5 3 1; 9 7 5 3 1; 13 7 5 3 1;
15 8 5 3 1; 14 9 5 3 1; 13 10 5 3 1; 14 13 5 3 1;
13 10 6 3 1; 14 10 6 3 1; 15 10 6 3 1; 11 9 7 3 1;
12 9 7 3 1; 13 11 7 3 1; 14 11 7 3 1; 15 11 7 3 1;
14 12 7 3 1; 14 11 10 3 1; 14 12 10 3 1; 14 14 11 3 1;
15 14 11 3 1; 15 10 6 4 1; 14 11 6 4 1; 13 12 7 4 1;
14 12 7 4 1; 15 12 7 4 1; 15 12 10 4 1; 14 12 11 4 1;
15 8 6 5 1; 14 11 8 5 1; 15 11 8 5 1; 14 12 8 5 1;
15 12 8 5 1; 13 10 8 6 1; 15 10 8 6 1; 14 11 8 6 1;
11 9 8 7 1; 13 11 8 7 1; 14 12 8 7 1; 14 12 10 8 1;
15 12 10 8 1; 14 11 11 8 1; 14 12 11 8 1; 14 14 11 8 1
];

tuple = dc4_tuple_gf32;

for row = 1:M
    col_id = find(H(row, :) ~= 0);
    tuple_index =randi([1,size(tuple,1)],1,1);
    [~, tupel_permutation] = sort(randi([1,1000000],size(col_id)));
    tupel_permutation = mod(tupel_permutation, size(tuple, 2));
    tupel_permutation(find(tupel_permutation == 0)) = randi([1,size(tuple, 2)], 1,sum(tupel_permutation == 0));
    tupel_multiplication = randi(qAry-1, 1, 1);
    vaule_rand = MULTIPLY_TABLE(tupel_multiplication + 1, tuple(tuple_index, tupel_permutation) + 1);
    H_rand(row, col_id) = vaule_rand;
end
% H_rand = randi(qAry-1, M,N);
% H = H .* H_rand;

%% write the nonbinary ldpc config into file
fid = fopen(NonBinaryFileName,'w');
fprintf(fid,'%d %d %d\n',N, M, qAry);
% fprintf(fid, 'not used\n');
% fprintf(fid, 'not used\n');
% fprintf(fid, 'not used\n');
fprintf(fid, '%d ', maxDegree);
fprintf(fid, '\n');
fprintf(fid, '%d ', colWei);
fprintf(fid, '\n');
fprintf(fid, '%d ', rowWei);
fprintf(fid, '\n');

for j=1:N
    for i=1:M
        if H(i,j) ~= 0
            fprintf(fid,'%d %d ',i, H_rand(i,j));
        end
    end
    fprintf(fid,'\n');
end

for i=1:M
    for j=1:N
        if H(i,j) ~= 0
            fprintf(fid,'%d %d ',j, H_rand(i,j));
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);