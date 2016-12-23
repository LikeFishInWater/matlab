    clear;
    x = logical([1 0 1 1 0 1 0 1 1 1 0 1]');
    hGen = comm.CRCGenerator([24 23 18 17 14 11 10 7 6 5 4 3 1 0]);
    codeword = step(hGen, x);
    y=CRCEncode(x,length(x),24,0);
    codeword-y
    hDetect = comm.CRCDetector([24 23 18 17 14 11 10 7 6 5 4 3 1 0]);
    [tx, err] = step(hDetect, codeword); 
    tx-x