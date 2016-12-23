function [H,rows,cols] = Gallager_gen_LDPC(miu,wr,wc)
% Code for generating regular parity-check matrix H by miu,wc,wr,using Gallager's approach
% 通过参数miu,wc,wr生成规则LDPC校验矩阵H，使用Gallager的办法
% Copyright sTeven joNes
% Email:sun.noon@gmail.com
% 2007-10-21

% Example:
% miu = 3;
% wr = 4;
% wc = 3;
% [H,rows,cols] = Gallager_gen_LDPC(miu,4608,,wc)

% example's OUTPUT:
% H =
%      1     1     1     1     0     0     0     0     0     0     0     0
% H(1) 0     0     0     0     1     1     1     1     0     0     0     0
%      0     0     0     0     0     0     0     0     1     1     1     1

%      0     1     0     0     1     0     0     0     0     1     1     0
% H(2) 0     0     1     1     0     1     1     0     0     0     0     0
%      1     0     0     0     0     0     0     1     1     0     0     1

%      0     0     0     1     0     1     1     1     0     0     0     0
% H(3) 0     0     1     0     1     0     0     0     1     1     0     0
%      1     1     0     0     0     0     0     0     0     0     1     1
% rows = 9
% cols = 12
% end of example

% H's size will be:[miu*wc,miu*wr]
% number of row    : miu*wc
% number of colume : miu*wr
% 矩阵H的尺寸为：[miu*wc,miu*wr]
% 行数：miu*wc
% 列数：miu*wr

% INPUT:
% miu : int（生成参数）
% wc  : weight of column（列重）（每列中1的个数）
% wr  : weight of row（行重10）（每行中1的个数）
% OUTPUT:
% H : Low-density parity-check matrix
% rows : 行数：miu*wc
% cols : 列数：miu*wr


%reference:
%[1]R.Gallager,"Low-density parity-check codes",MIT press,1963
%[2]Willam E. Ryan,"An Introduction to LDPC Codes"

rows = miu*wc;
cols = miu * wr;
h = zeros(rows,cols);

wrones = ones(1,wr);

% H = [H(1)
%      H(2)
%      .
%      .
%      .
%      H(wc)] 

for ii = 1 : miu                          %完成第一个矩阵的构造
    h(ii,(ii-1)*wr + 1:ii*wr) = wrones;
end
% H(1)...done

%算法：对于每一行，随机选取第i列，如果这一列的标志位（flag(i)）为0，
%说明在该子矩阵H(topi)中，可以为该行的第i列赋值为1，当一行中赋值达到行重wr时，开始对下一行进行处理
for topi = 2:wc%for H(2) to H(wc) 对第二个到第wc个矩阵
    flag  = zeros(1,cols);%该H(topi)中的每一行是否已经分配了1的标志
    for ii = (topi-1)*miu +1 : topi * miu % for each line of H(topi)
        countones_eachline = 0;%对于每一行中的1的个数，进行计数，当1的个数等于行重wr时，进行下一行的处理
        while (countones_eachline<wr)
            index_one_col = ( round(rand(1) * (cols-1)) )+1;%随机选取一列rand1产生01之间的均匀数。  1+0.11            
            if (flag(index_one_col) ~= 1)%标志位还为0的话，该列还未设置过，可以进行下去
                h(ii,index_one_col) = 1;%矩阵中的相应选取到的元素设为1
                flag(index_one_col) = 1;%标志位置1
                countones_eachline = countones_eachline + 1;
            end%end if (flag(index_one_col) ~= 1)             
        end%end while (countones_eachline<=wr)
    end%end for ii = (topi-1)*miu +1 : topi * miu  
end%end for topi = 2:wc
t=1;
for i=1:miu*wr
    for j=1:miu*wc
        if h(j,i)==1
            H_col_id(t)=i;
            H_row_id(t)=j;
            t=t+1;
        end
    end 
end
% H = h;
H=sparse(H_row_id,H_col_id,ones(1,length(H_row_id)));

% you can get H(n) in your code as below
% H1 = h(1:miu,1:cols);
% Hn = h((n-1)*miu +1:n * miu,1:cols)