function rateMatchOut = rate_match( rateMatchLen,codeBlockSize,codeBlockIndex,nullBitNum,systemBit,parityBit,iParityBit,R,C,k0 )
%UNTITLED5 此处显示有关此函数的摘要
%   此处显示详细说明
    %%%% turbo-encoder output %%%%
    codeBlockSizeAddTail=codeBlockSize+4;

    %%%% sub-block interleaver %%%%
    interMatrixSize=R*C;
    cycleBufSize=3*interMatrixSize;
    F0=interMatrixSize-codeBlockSizeAddTail;
    P=[0,16,8,24,4,20,12,28,2,18,10,26,6,22,14,30,1,17,9,25,5,21,13,29,3,19,11,27,7,23,15,31]+1;
    %P=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32];

    %%%% sub-block interleaver and bit collection starting %%%%
    %input :3 channel
    %output:cycleBuf
    temp1 = sub_interleave(systemBit,0,R,C,F0);
    cycleBufSaveAddr=1;
    cycleBuf=zeros(1,3*R*C);
    for i=1:R*C
        cycleBuf(cycleBufSaveAddr) = temp1(i);
        cycleBufSaveAddr=cycleBufSaveAddr+1;
    end
    temp2 = sub_interleave(parityBit,0,R,C,F0);
    for i=1:R*C
        cycleBuf(cycleBufSaveAddr) = temp2(i);
        cycleBufSaveAddr=cycleBufSaveAddr+2;
    end
    temp3 = sub_interleave(iParityBit,1,R,C,F0);
    cycleBufSaveAddr = R*C + 2;
    for i=1:R*C
        cycleBuf(cycleBufSaveAddr) = temp3(i);
        cycleBufSaveAddr=cycleBufSaveAddr+2;
    end
    %%%% sub-block interleaver and bit collection ending %%%%

    %%%% rate match starting %%%%
    cycleBufGetAddr = k0;%set initial get address
    rateMatchSaveAddr=1;
    rateMatchOut=zeros(1,rateMatchLen);
    while (rateMatchSaveAddr <= rateMatchLen)
        % decide which bit stream it is , 0->system bit, 1->parity bit, 2-> Turbo_Interleave parity bit
        if (cycleBufGetAddr <= interMatrixSize)
            bitType = 0;
            bitPosition = cycleBufGetAddr;
        else
            if mod((cycleBufGetAddr-interMatrixSize),2) == 0
                bitType=2;
            else 
                bitType=1;
            end
            bitPosition = ceil((cycleBufGetAddr-interMatrixSize)/2);
        end
        % add another null bit
        if (0 == codeBlockIndex && 2 ~= bitType)
            totalNullBitNum = F0 + nullBitNum;
        else
            totalNullBitNum = F0;
        end
        if bitType==2
            offset=1;
        else 
            offset=0;
        end
        originBitPosition = mod((P(ceil(bitPosition/R)) + C*(mod(bitPosition-1,R))+offset)-1,interMatrixSize)+1;         
        if (originBitPosition > totalNullBitNum)
            rateMatchOut(rateMatchSaveAddr) = cycleBuf(cycleBufGetAddr);
            rateMatchSaveAddr=rateMatchSaveAddr+1;
        end
        cycleBufGetAddr = mod((cycleBufGetAddr+1),cycleBufSize);
    end
    %%%% rate match ending %%%%


end

