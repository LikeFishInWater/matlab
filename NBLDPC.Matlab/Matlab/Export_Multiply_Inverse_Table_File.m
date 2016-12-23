clear; clc;

%% File Name
TableFileName = './matrix/Table.GF.16.txt';
qAry = 16;

global MULTIPLY_TABLE;
global INVERSE_TABLE;

%% Generate the GF(q) table
NonBinary_LDPC_GenerateTable(qAry);

%% Export the Table to File
fid = fopen(TableFileName,'w');
fprintf(fid,'Multiply Table:\n');
for j=1:qAry
    for i=1:qAry
        fprintf(fid,'%d ',MULTIPLY_TABLE(i,j));
    end
    fprintf(fid,'\n');
end

fprintf(fid,'Inverse Table:\n');
for i=1:qAry-1
    fprintf(fid,'%d ',INVERSE_TABLE(i));
end
fprintf(fid,'\n');
fclose(fid);