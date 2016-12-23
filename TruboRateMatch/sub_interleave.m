function  outputBit = sub_interleave( inputBit,offset,R,C,F0 )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
%   C=32;
    P=[0,16,8,24,4,20,12,28,2,18,10,26,6,22,14,30,1,17,9,25,5,21,13,29,3,19,11,27,7,23,15,31]+1;
    %P=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
    getAddr=1;
    saveAddr=1;
    rowIndex=1; 
    columnIndex=1;
    
    while saveAddr <= R*C
        getAddr = mod((C*(rowIndex-1) + P(columnIndex) + offset)-1,(R*C))+1;
        if (getAddr > F0)
			outputBit(saveAddr) = inputBit(getAddr-F0);
        else
			outputBit(saveAddr) = 0; %%%fill null bit
        end
		saveAddr=saveAddr+1;
		if (rowIndex < R)
			rowIndex=rowIndex+1;
        else
			rowIndex = 1;
			columnIndex=columnIndex+1;
        end
    end
end

