%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NAME:        NonBinary_LDPC_Encoder_Check
% PURPOSE:     Check the whether the encoded code is right
%
% Input:
% H: GF(qAry) NonBinary LDPC check Matrix
% R: R is the indexes of variable nodes connted to Chk_j
% C: C is the indexes of check nodes connted to Var_i
% qAry: GF(qAry)
% EncodeInput: the random msg sequence, which is uncoded
% EncodeOutput: the coded sequence, which is coded
%
% Output:
% isCorrect: 1:right 0:error
%
% AUTHOR:       Xiaoshi
% DATE:         2014.12.10
% VERSION:      v1.0
% REVISED BY:   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function isCorrect = NonBinary_LDPC_Encoder_Check(H, R, C, qAry, EncodeInput, EncodeOutput)
[check_size, var_size] = size(H);

isCorrect = true;
syndrome = zeros(check_size, 1);

if EncodeOutput(1,1:var_size - check_size) ~= EncodeInput
    isCorrect = false;
else
    for j=1:check_size
        [~, link_sum] = size(R{j});
        syndrome(j, 1) = 0;
        for var_no=1:link_sum
            i = R{j}(var_no);
            syndrome(j,1) = nonbinary_add(syndrome(j,1), nonbinary_multiply(EncodeOutput(1, i), H(j,i))); 
        end
    end
    if sum(syndrome(:,1)) ~= 0 
        isCorrect = false;
    end
end