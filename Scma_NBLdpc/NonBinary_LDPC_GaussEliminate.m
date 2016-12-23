%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NAME:        NonBinary_LDPC_GaussEliminate
% PURPOSE:     Gauss Eliminate Process of matrix H
%
% Input:
% H: the original matrix H
% qAry: GF(qAry)
%
% Output:
% Eliminated_H:  Eliminated H
% R: the indexes of variable nodes connted to Chk_j in Eliminated_H
% C: the indexes of check nodes connted to Var_i in Eliminated_H
% ExchangedOriginCol: the original col index need to be exchanged
% ExchangedDestinCol: the destination col index need to be exchanged
%
% AUTHOR:       Xiaoshi
% DATE:         2014.12.10
% VERSION:      v1.0
% REVISED BY:   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [Eliminated_H, R, C, ExchangedOriginCol, ExchangedDestinCol] = NonBinary_LDPC_GaussEliminate(H, qAry)
[M, N] = size(H);
ExchangedOriginCol =[];
ExchangedDestinCol=[];
for row = M:-1:1
    col = N - M + row;
    %if H(row, col) = 0, which is the last element of this row, then seek
    %from bottom to the top, until H(seek_row, col) is an nonbinary element
    if H(row, col) == 0
        isfound = false;
        for seek_row = row - 1:-1:1
            if H(seek_row, col) ~= 0
                isfound = true;
                tempExchange = H(row, :);
                H(row, :) = H(seek_row, :);
                H(seek_row, :) = tempExchange;
                break;
            end
        end
        if isfound == false
            for seek_col = col - 1:-1:1
                if H(row, seek_col) ~= 0
                    isfound = true;
                    ExchangedOriginCol = [col ExchangedOriginCol];
                    ExchangedDestinCol = [seek_col ExchangedDestinCol];
                    tempExchange = H(:, seek_col);
                    H(:, seek_col) = H(:, col);
                    H(:, col) = tempExchange;
                    break;
                end
            end
        end
    end
    for seek_row = row - 1:-1:1
        if H(seek_row, col) ~= 0
            k = nonbinary_multiply(nonbinary_inverse(H(seek_row, col)), H(row, col));
            for var_index = 1:col
                H(seek_row, var_index) = nonbinary_add(nonbinary_multiply(H(seek_row, var_index), k), full(H(row, var_index))); 
            end
        end       
    end
end
Eliminated_H = H;

%Generate the link version
R=cell(M,1);
C=cell(1,N);
for j=1:M
    R{j} = find(Eliminated_H(j,:)~=0);
end

for i=1:N
    C{i} = find(Eliminated_H(:,i)~=0);
end