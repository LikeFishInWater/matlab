function [ systemBit,parityBit,iParityBit ] = tiqu( encodedData)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
    dataLen=length(encodedData);
    for i=1:dataLen/3
        systemBit(i)=encodedData(3*i-2);
        parityBit(i)=encodedData(3*i-1);
        iParityBit(i)=encodedData(3*i);
    end
end

