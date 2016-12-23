%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NAME:        nonbinary_add
% PURPOSE:     NonBinary LDPC decoder ADD function module qAry
%
% Input:
% input1: GF(q) element
% input2: GF(q) element
%
% Output:
% output: output = input1 + input2, this operator is on GF(q) field
%
% AUTHOR:       Xiaoshi
% DATE:         2014.12.10
% VERSION:      v1.0
% REVISED BY:   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function output = nonbinary_add(input1, input2)
output = bitxor(input1, input2);