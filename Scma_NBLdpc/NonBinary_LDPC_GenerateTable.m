%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NAME:        NonBinary_LDPC_GenerateTable
% PURPOSE:     Generate the */^{-1} table of GF(qAry)
%
% Input:
% qAry: GF(qAry)
%
% AUTHOR:       Xiaoshi
% DATE:         2014.12.04
% VERSION:      v1.0
% REVISED BY:   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function NonBinary_LDPC_GenerateTable(qAry)
global MULTIPLY_TABLE;
global INVERSE_TABLE;

MULTIPLY_TABLE = zeros(qAry, qAry);
INVERSE_TABLE = zeros(qAry - 1, 1);

for row_a = 1:qAry
    for col_b = 1:qAry
        gf_a = gf(row_a - 1, log2(qAry));
        gf_b = gf(col_b - 1, log2(qAry));
        for q=0:qAry-1
            gf_q = gf(q, log2(qAry)); 
            if gf_a * gf_b == gf_q
                MULTIPLY_TABLE(row_a, col_b) = q;
                break;
            end
        end
    end
end

for row_a = 2:qAry
    for col_b = 2:qAry
        if MULTIPLY_TABLE(row_a, col_b) == 1
            INVERSE_TABLE(row_a - 1) = col_b - 1;
            break;
        end
    end
end

% switch qAry
%     case 4,
%         MULTIPLY_TABLE = [0 0 0 0; 0 1 2 3; 0 2 3 1; 0 3 1 2];
%         INVERSE_TABLE = [1;3;2];
%     case 8,
%         MULTIPLY_TABLE = [0 0 0 0 0 0 0 0; 0 1 2 3 4 5 6 7; 0 2 4 6 3 1 7 5; 0 3 6 5 7 4 1 2;
%             0 4 3 7 6 2 5 1; 0 5 1 4 2 7 3 6; 0 6 7 1 5 3 2 4; 0 7 5 2 1 6 4 3];
%         INVERSE_TABLE = [1;5;6;7;2;3;4];
%     case 16,
%         MULTIPLY_TABLE = [0 0 0 0; 0 1 2 3; 0 2 3 1; 0 3 1 2];
%         INVERSE_TABLE = [1;3;2];
%     case 32,
%         MULTIPLY_TABLE = [0 0 0 0; 0 1 2 3; 0 2 3 1; 0 3 1 2];
%         INVERSE_TABLE = [1;3;2];
%     case 64,
%         MULTIPLY_TABLE = [0 0 0 0; 0 1 2 3; 0 2 3 1; 0 3 1 2];
%         INVERSE_TABLE = [1;3;2];
%     case 128,
%         MULTIPLY_TABLE = [0 0 0 0; 0 1 2 3; 0 2 3 1; 0 3 1 2];
%         INVERSE_TABLE = [1;3;2];
%     case 256
%         MULTIPLY_TABLE = [0 0 0 0; 0 1 2 3; 0 2 3 1; 0 3 1 2];
%         INVERSE_TABLE = [1;3;2];
%         
% end
