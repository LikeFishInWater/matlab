%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NAME:        nonbinary_multiply
% PURPOSE:     NonBinary LDPC decoder multiply function input1 * input2
%
% Input:
% input1: GF(q) element
% input2: GF(q) element
%
% Output:
% output: output = input1 * input2, this operator is on GF(q) field
%
% AUTHOR:       Xiaoshi
% DATE:         2014.12.10
% VERSION:      v1.0
% REVISED BY:   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function output = nonbinary_multiply(input1, input2)
global MULTIPLY_TABLE;
output = MULTIPLY_TABLE(input1 + 1, input2 + 1);