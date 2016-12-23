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
function NonBinary_LDPC_ConvertMatrixFile(BinaryFileName, NonBinaryFileName, qAry)
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
H = (H >= 1);
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
H_rand = randi(qAry-1, M,N);
H = H .* H_rand;

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