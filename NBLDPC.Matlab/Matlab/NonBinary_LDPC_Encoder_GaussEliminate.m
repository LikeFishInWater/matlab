%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NAME:        Nonbinary LDPC Encode
% PURPOSE:     Eliminated_H is a Gauss Eliminated maxtrix
%
% Input:
% EncodeInput: 1*(N-M)  Msg Vector 
% Eliminated_H: Eliminated GF(qAry) NonBinary LDPC check Matrix
% GaussEliminate_R: the indexes of variable nodes connted to Chk_j
% GaussEliminate_C: the indexes of check nodes connted to Var_i
% ExchangedOriginCol: the original col index need to be exchanged
% ExchangedDestinCol: the destination col index need to be exchanged
% qAry: GF(qAry)
%
% Output:
% EncodeOutput: 1*N Encoded Vector 
%
% AUTHOR:       Xiaoshi
% DATE:         2014.12.10
% VERSION:      v1.0
% REVISED BY:   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function EncodeOutput = NonBinary_LDPC_Encoder_GaussEliminate(EncodeInput, Eliminated_H, GaussEliminate_R, GaussEliminate_C, ExchangedOriginCol, ExchangedDestinCol, qAry)
[M,N] = size(Eliminated_H);
EncodeOutput = zeros(1,N);
EncodeOutput(1,1:N-M) = EncodeInput;

for j=1:M
    col = j + N - M;
    [~,link_sum] = size(GaussEliminate_R{j});
    if GaussEliminate_R{j}(link_sum) ~= col
        continue;
    end
    for var_no = 1:link_sum-1
        i = GaussEliminate_R{j}(var_no);
        EncodeOutput(1,col) = nonbinary_add(EncodeOutput(1,col), nonbinary_multiply(EncodeOutput(1,i), Eliminated_H(j,i)));
    end
    EncodeOutput(1,col) = nonbinary_multiply(EncodeOutput(1,col), nonbinary_inverse(Eliminated_H(j,col)));
end

% for j = 1:M
%     col = j + N - M;
%     for i=1:col-1
%         EncodeOutput(1,col) = nonbinary_add(EncodeOutput(1,col), nonbinary_multiply(EncodeOutput(1,i), Eliminated_H(j,i)));
%     end
%     EncodeOutput(1,col) = nonbinary_multiply(EncodeOutput(1,col), nonbinary_inverse(Eliminated_H(j,col)));
% end

[~,exchange_size] = size(ExchangedOriginCol);
for i = 1:exchange_size
    temp = EncodeOutput(1, ExchangedOriginCol(i));
    EncodeOutput(1, ExchangedOriginCol(i)) = EncodeOutput(1, ExchangedDestinCol(i));
    EncodeOutput(1, ExchangedDestinCol(i)) = temp;
end