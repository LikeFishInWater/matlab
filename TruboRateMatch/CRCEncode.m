function seqOut = CRCEncode(seqIn,seqInLen,crcLen,crc24Type)
	% -	gCRC24A(D) = [D24 + D23 + D18 + D17 + D14 + D11 + D10 + D7 + D6 + D5 + D4 + D3 + D + 1] and;
	% -	gCRC24B(D) = [D24 + D23 + D6 + D5 + D + 1] for a CRC length L = 24 and;
	% -	gCRC16(D) = [D16 + D12 + D5 + 1] for a CRC length L = 16.
	% -	gCRC8(D) = [D8 + D7 + D4 + D3 + D + 1] for a CRC length of L = 8.
    
%     clear;
%     x = logical([1 0 1 1 0 1 0 1 1 1 0 1]');
%     hGen = comm.CRCGenerator([8 7 4 3 1 0]);
%     codeword = step(hGen, x);
%     y=CRCEncode(x,length(x),8,0);
%     codeword-y
    reg=zeros(1,24);
    G=zeros(1,25);
    switch(crcLen)
        case 8,
            G(0+1)=1;
            G(1+1)=1;
            G(3+1)=1;
            G(4+1)=1;
            G(7+1)=1;
            G(8+1)=1;
        case 16,
            G(0+1)=1;
            G(5+1)=1;
            G(12+1)=1;
            G(16+1)=1;
        case 24,
            if (crc24Type==0)
                G=[1,1,0,1,1,1,1,1,0,0,1,1,0,0,1,0,0,1,1,0,0,0,0,1,1];
            else
                G=[1,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1];
            end
    end
    if (crcLen ~= 0)
		for i=1:seqInLen %(i = 0; i < seqInLen; i++)
			temp = xor(reg(crcLen),seqIn(i));
			for j=crcLen:-1:2
				reg(j) = xor(reg(j -1) , (G(j) && temp));
            end
			reg(1) = temp;
        end
        seqOut=seqIn;
		for i=1:crcLen
			seqOut(i+seqInLen) = reg(crcLen+1-i); 
        end
    else  % if CRC length is 0, just copy input to output
		if (seqIn ~= seqOut)
            seqOut=seqIn;
        end
    end
end