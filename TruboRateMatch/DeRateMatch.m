function [ systemBit,parityBit,iParityBit ] = DeRateMatch( decodeIn,blockSizeAddTail,k0,R,C,rateMatchLen,nullBitNum,codeBlockIndex )
%UNTITLED6 此处显示有关此函数的摘要
%   此处显示详细说明
    P=[0,16,8,24,4,20,12,28,2,18,10,26,6,22,14,30,1,17,9,25,5,21,13,29,3,19,11,27,7,23,15,31]+1;
    %P=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
    interMatrixSize = R * C;
    F0 = interMatrixSize - blockSizeAddTail;
    getAddr=k0;
    getBitLen=1;
    maxInfoLength=blockSizeAddTail-4;
    systemBit=zeros(1,maxInfoLength+4);
    parityBit=zeros(1,maxInfoLength+4);
    iParityBit=zeros(1,maxInfoLength+4);
    sysLen=0;
    parLen=0;
    iParLen=0;
    while (getBitLen <= rateMatchLen)
		%decide which bit stream it is , 0->system bit, 1->parity bit, 2-> interleave parity bit
		if (getAddr <= interMatrixSize)
			bitType = 0;
			bitPosition = getAddr; %系统bits
        else
            if mod((getAddr-interMatrixSize),2) == 1
                bitType=1;
            else
                bitType=2;
            end
			bitPosition = ceil((getAddr-interMatrixSize)/2);%校验bits
        end
        if bitType==2
            offset=1;
        else 
            offset=0;
        end
		originBitPosition = mod((P(ceil(bitPosition/R)) + C*(mod(bitPosition-1,R)) + offset)-1,interMatrixSize)+1;
		if (0 == codeBlockIndex && 2 ~= bitType)
			totalNullBitNum = F0 + nullBitNum;  %for sub-block 0's system and parity bit stream, add additional nullBitNum NULL bits.
        else
			totalNullBitNum = F0;
        end

		if (originBitPosition > F0)
			switch (bitType)
			case 0,
				if (originBitPosition > totalNullBitNum )
					systemBit(originBitPosition-F0) = decodeIn(getBitLen);
					getBitLen=getBitLen+1;
					sysLen=sysLen+1;
                else
					systemBit(originBitPosition-F0) = 0;%null bit added to first code block
                end
			case 1,
				if (originBitPosition > totalNullBitNum )
					parityBit(originBitPosition-F0) = decodeIn(getBitLen);
					getBitLen=getBitLen+1;
					parLen=parLen+1;
                else
					parityBit(originBitPosition-F0) = 0;%null bit added to first code block
                end
			case 2,
				iParityBit(originBitPosition-F0) = decodeIn(getBitLen);
				getBitLen=getBitLen+1;
				iParLen=iParLen+1;
            end
        end
		getAddr =mod((getAddr+1),(3*interMatrixSize));
    end


end

