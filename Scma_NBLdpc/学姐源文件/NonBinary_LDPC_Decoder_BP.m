%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NAME:        Nonbinary LDPC Decode
% PURPOSE:     Universal LLR BP decoder for Nonbinary LDPC
%
% Input:
% H: GF(qAry) NonBinary LDPC check Matrix
% R: R is the indexes of variable nodes connted to Chk_j
% C: C is the indexes of check nodes connted to Var_i
% qAry: GF(qAry)
% L_ch: (qAry-1) * N matrix, L_ch(:,i) represents the LLR of variable i
% max_iterations: the maximum iteration time
%
% Output:
% DecodeOutput: 1*N Decoded Vector 
%
% AUTHOR:       Xiaoshi
% DATE:         2014.12.04
% VERSION:      v1.1
% REVISED BY:   Xiaoshi--2015.01.27--Revised the LLRBoxPlus function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function DecodeOutput = NonBinary_LDPC_Decoder_BP(H, R, C, qAry, L_ch, max_iterations)
%C,Column每列所连校验节点的编号
[check_size, var_size] = size(H);
iterations=0; 
%% initialization
q = zeros(qAry-1, check_size, var_size);
r = zeros(qAry-1, check_size, var_size);
for i = 1:var_size
    [link_sum, ~] = size(C{i});
    for check_no = 1:link_sum
        j = C{i}(check_no);
        q(:,j,i) = L_ch(:,i);
    end
end
for j = 1:check_size
    [~, link_sum] = size(R{j});
    for var_no = 1:link_sum
        i = R{j}(var_no);
        r(:,j,i) = zeros(qAry-1, 1);
    end
end

%% iterative decoding process
L_post = zeros(qAry-1, var_size);
DecodeOutput = zeros(1, var_size);
syndrome = zeros(check_size, 1);
while true   
    %% tentative decoding
    for i=1:var_size
        L_post(:,i) = L_ch(:,i);
        [link_sum, ~] = size(C{i});
        for check_no = 1:link_sum
            j = C{i}(check_no);
            L_post(:,i) = L_post(:,i) + r(:,j,i);
        end
        [maxLLR, alpha_i]= max(L_post(:,i));
        if maxLLR > 0
            DecodeOutput(1, i) = alpha_i;
        else
            DecodeOutput(1, i) = 0;
        end
    end
    %% check parity 
    for j=1:check_size
        [~, link_sum] = size(R{j});
        syndrome(j, 1) = 0;
        for var_no=1:link_sum
            i = R{j}(var_no);
            syndrome(j,1) = nonbinary_add(syndrome(j,1), nonbinary_multiply(DecodeOutput(1, i), H(j,i))); 
        end
    end
    if sum(syndrome(:,1)) == 0 || iterations >= max_iterations
        break;
    else
        iterations = iterations + 1;
    end
    %% Horizontal step
    for i=1:var_size
        [link_sum, ~] = size(C{i});
        log_sum = 0;
        for check_no = 1:link_sum
            j = C{i}(check_no);
            log_sum = log_sum + r(:,j,i);
        end
%         log_sum = sum(r(:,:,i),2);
        for check_no=1:link_sum
            j = C{i}(check_no);
            q(:,j,i) = L_ch(:,i) + log_sum - r(:,j,i);
        end
    end
    %% Vertical step
    for j=1:check_size
        [~, link_sum]  = size(R{j});
        for var_no = 1:link_sum
            %forward
            L_sigma = zeros(qAry-1,1);
            for l = 1 : var_no - 1
                i = R{j}(l);
                if l == 1
                    L_sigma = LLR_BoxPlus(L_sigma, q(:,j,i), 0, H(j,i), qAry);
                else
                    L_sigma = LLR_BoxPlus(L_sigma, q(:,j,i), 1, H(j,i), qAry);
                end
            end
            %backward
            L_rho = zeros(qAry-1,1);
            for l = link_sum : -1 : var_no + 1
                i = R{j}(l);
                if l == link_sum
                    L_rho = LLR_BoxPlus(L_rho, q(:,j,i), 0, H(j,i), qAry);
                else
                    L_rho = LLR_BoxPlus(L_rho, q(:,j,i), 1, H(j,i), qAry);
                end
            end
            %given message
            if var_no == 1
                i = R{j}(var_no);
                r(:,j,i) = LLR_BoxPlus(L_sigma, L_rho, 0, nonbinary_inverse(H(j,i)), qAry);
            elseif var_no == link_sum
                i = R{j}(var_no);
                r(:,j,i) = LLR_BoxPlus(L_sigma, L_rho, nonbinary_inverse(H(j,i)), 0, qAry);
            else
                i = R{j}(var_no);
                r(:,j,i) = LLR_BoxPlus(L_sigma, L_rho, nonbinary_inverse(H(j,i)), nonbinary_inverse(H(j,i)), qAry);
            end
        end
    end
end

%% NonBinary LDPC decoder Box Plus function
function L = LLR_BoxPlus(L1, L2, A1, A2, qAry)
if A1 == 0
    for alpha_i = 1:qAry - 1
        L(alpha_i,1) = L2(nonbinary_multiply(nonbinary_inverse(A2), alpha_i));
    end
elseif A2 == 0
    for alpha_i = 1:qAry - 1
        L(alpha_i,1) = L1(nonbinary_multiply(nonbinary_inverse(A1), alpha_i));
    end
else
    for alpha_i = 1:qAry - 1
        a = exp(L1(nonbinary_multiply(nonbinary_inverse(A1), alpha_i)));
        b = exp(L2(nonbinary_multiply(nonbinary_inverse(A2), alpha_i)));
        c = 0;
        for x = 1:qAry-1
            if x ~= nonbinary_multiply(alpha_i, nonbinary_inverse(A1))
                c = c + exp(L1(x) + L2(nonbinary_multiply(nonbinary_inverse(A2) , nonbinary_add(alpha_i, nonbinary_multiply(x, A1)))));
            end
        end
        d = 0;
        for x = 1:qAry-1
            d = d + exp(L1(x) + L2(nonbinary_multiply(nonbinary_inverse(A2) , nonbinary_multiply(x, A1))));
        end
        L(alpha_i,1) = log(a + b + c) - log(1 + d);
    end
end


% for alpha_i = 1:qAry - 1
%     a = L1(nonbinary_multiply(nonbinary_inverse(A1), alpha_i));
%     b = L2(nonbinary_multiply(nonbinary_inverse(A2), alpha_i));
%     c = max(a,b);
%     for x = 1:qAry-1
%         if x ~= nonbinary_multiply(alpha_i, nonbinary_inverse(A1))
%             c = max(c, L1(x) + L2(nonbinary_multiply(nonbinary_inverse(A2) , nonbinary_add(alpha_i, nonbinary_multiply(x, A1)))));
%         end
%     end
%     d = 0;
%     for x = 1:qAry-1
%         d = max(d, L1(x) + L2(nonbinary_multiply(nonbinary_inverse(A2) , nonbinary_multiply(x, A1))));
%     end
%     L(alpha_i,1) = c - d;
% end
