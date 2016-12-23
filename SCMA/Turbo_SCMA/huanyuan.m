function [ encodedData ] = huanyuan( systemBit,parityBit,iParityBit )
%UNTITLED7 此处显示有关此函数的摘要
%   此处显示详细说明
    for i=1:length(systemBit)
        encodedData(3*i-2)=systemBit(i);
        encodedData(3*i-1)=parityBit(i);
        encodedData(3*i)=iParityBit(i);
    end
end

