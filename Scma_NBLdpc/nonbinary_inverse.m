%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NAME:        nonbinary_inverse
% PURPOSE:     NonBinary LDPC decoder inverse function input^-1 = 1 / input
%
% Input:
% input: GF(q) element
%
% Output:
% output: output = input^-1 = 1 / input, this operator is on GF(q) field
%
% AUTHOR:       Xiaoshi
% DATE:         2014.12.10
% VERSION:      v1.0
% REVISED BY:   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function output = nonbinary_inverse(input)
global INVERSE_TABLE
if input ~= 0
    output = INVERSE_TABLE(input);
else
    disp('Div 0 Error');
end