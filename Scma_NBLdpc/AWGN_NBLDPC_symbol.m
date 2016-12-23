clear;
clc;
%% Initial Simulation Parameters
MaxIter =50;
nsample =100000;
[qAry_cpp, TABLE_MULTIPLY_cpp, TABLE_ADD_cpp, TABLE_INVERSE_cpp, maxVarDegree_cpp, maxChkDegree_cpp, VarDegree_cpp, ChkDegree_cpp, ...
    H_cpp, Row_Link_Col_cpp, Col_Link_Row_cpp,Q_cpp, R_cpp, L_Post_cpp, LDR_Vector_cpp, L_SIGMA_cpp, L_RHO_cpp,CodeLen_cpp, ChkLen_cpp, MsgLen_cpp, ...
    Rate_cpp] =Initial_CPP('./matrix/BM.N28.M16.Z139.C50.G10.tcj.1.scale.4.tanner.GF.8.txt', './matrix/Arith.Table.GF.8.txt'); %Initial_CPP('./matrix/PEG.1334.GF.64.txt', './matrix/TABLE.GF64.txt');

%% Initial NBLDPC
NonBinaryFileName ='./matrix/BM.N28.M16.Z139.C50.G10.tcj.1.scale.4.tanner.GF.8.txt';
global qAry
qAry =8;
NonBinary_LDPC_GenerateTable(qAry);
[H, N, M, Link_Var_to_Chk, Link_Chk_to_Var, qAry] = NonBinary_LDPC_LoadMatrixFile(NonBinaryFileName);
[Eliminated_H, GaussEliminate_R, GaussEliminate_C, ExchangedOriginCol, ExchangedDestinCol]  = NonBinary_LDPC_GaussEliminate(H, qAry);
R = (N - M) / N;

%% Initial Modulator & Demodulator
global ModType
ModType=4;
hMod= comm.RectangularQAMModulator('ModulationOrder',2^ModType,'BitInput',true,'NormalizationMethod',...
    'Average power');

%% Initial MIMO Parameter
Eb_N0_dB=[2:0.5:7];
Eb_N0= 10.^(Eb_N0_dB./10);
noise_real_sigma2_vec = 0.5/ModType/R./Eb_N0;

%% Initial Others
BER=zeros(1,length(Eb_N0_dB)) ;
FER=zeros(1,length(Eb_N0_dB)) ;
err_bit=zeros(nsample,length(Eb_N0_dB));
err_bit_sum=zeros(1,length(Eb_N0_dB)) ;
err_Frames=zeros(1,length(Eb_N0_dB)) ;
countCycles=zeros(1,length(Eb_N0_dB)) ;
Symbol_ModType = ArrayModType(ModType, hMod);

q=log2(qAry);
min_lcm=lcm(q,ModType);
qam_piece=min_lcm/ModType;
qary_piece=min_lcm/q;
%if q>ModType
%    [need_bit,array_0]=qAry_to_QAM_numlist(q,ModType);
%end

%% Simulation
for Eb_N0_i=1:length(Eb_N0_dB)
    Eb_N0_dB(Eb_N0_i)
    noise_real_sigma2= noise_real_sigma2_vec (Eb_N0_i);
    for ns=1:nsample
        % NBLDPC Encoder      
        EncodeInput =randi(qAry, 1, N-M) - 1;
        EncodeOutput = NonBinary_LDPC_Encoder_GaussEliminate(EncodeInput, Eliminated_H, GaussEliminate_R, GaussEliminate_C, ExchangedOriginCol, ExchangedDestinCol, qAry);
        % Modulation
        %check_length(EncodeOutput,min_lcm);
        [Symbol,length_now] = LDPC2QAM(EncodeOutput,hMod, qAry,min_lcm);
        %Pr_simple=zeros(qAry,1);
        Pr_qary=zeros(qAry-1,length(Symbol));
        noise_sigma2=2*noise_real_sigma2;

        block=length_now/ModType;
        noise = wgn(block,1, 10 * log10(noise_sigma2), 'complex');
        Y= Symbol+ noise;
        for i = 1 : block/qam_piece
            Pr_simple_qam=zeros(2^ModType,qam_piece);
            for n=1:qam_piece
                for k = 1:2^ModType
                    Pr_simple_qam(k,n) = exp( -(abs(Y((i-1)*qam_piece+n) -  Symbol_ModType(k,1)).^2) / noise_sigma2);
                end
            end
            
            
            %%这一段是枚举法做出来的，不算通用程序
            switch q
                case 6
                    switch ModType
                        case 4
                            %need_bit=[4,2,4;0,2,0]

                            for l=0:qAry-1
                               ll=dec2bin(l,q);
                               a=bin2dec(ll(1:4))+1;
                               b= (bin2dec(ll(5:6))*4+1):  (bin2dec(ll(5:6))*4+4);
                               c=(bin2dec(ll(1:2))+1):4:(bin2dec(ll(5:6))+13);
                               d=bin2dec(ll(3:6))+1;
                               Pr_LLR_temp(l+1,1)=Pr_simple_qam(a,1)*sum(Pr_simple_qam(b,2) );
                               Pr_LLR_temp(l+1,2)=sum(Pr_simple_qam(c,2))*Pr_simple_qam(d,3);
                            end
                            
                            
                            
                            for u=1:qary_piece
                                for v=1:qAry-1
                                    Pr_LLR_two(v,u)=log(Pr_LLR_temp(v+1,u)/Pr_LLR_temp(1,u ));
                                end
                            end
                    end
                    e=(1+(i-1)*qary_piece):(i*qary_piece);
                            Pr_LLR(:,e)=Pr_LLR_two;
                            Pr_LLR_temp=zeros;
                            Pr_LLR_two=zeros;
                            
                case 8
                    switch ModType
                        case 4
                            need_bit=[4,4]
                        case 6
                            need_bit=[6,2,4,6;0,4,2,0]
                            
                    end
            end
            
            %%
            
           % Pr_LLR(:,(1+(i-1)*qary_piece):(i*qary_piece))=QAM_to_qAry(array_0,need_bit,Pr_simple_qam);
    
        end

        [DecodeOutput_uint, DecodeOutput]= Decode_CPP(Pr_LLR, CodeLen_cpp, ChkLen_cpp, VarDegree_cpp, ChkDegree_cpp, ...
            Col_Link_Row_cpp, Row_Link_Col_cpp, qAry_cpp, TABLE_MULTIPLY_cpp, TABLE_ADD_cpp, TABLE_INVERSE_cpp, ...
            H_cpp, Q_cpp, R_cpp, L_Post_cpp, LDR_Vector_cpp, L_SIGMA_cpp, L_RHO_cpp, MaxIter);
        
        err_bit(ns, Eb_N0_i) = biterr(EncodeOutput, DecodeOutput);
        if err_bit(ns, Eb_N0_i)  ~= 0
            err_Frames(Eb_N0_i) = err_Frames(Eb_N0_i) + 1;
            err_bit_sum(Eb_N0_i) = err_bit(ns,Eb_N0_i) + err_bit_sum(Eb_N0_i) ;
        end
        countCycles(Eb_N0_i) = countCycles(Eb_N0_i) + 1;
        disp(['Eb/N0: ' num2str(Eb_N0_dB(Eb_N0_i)) ' BER:' num2str(err_bit_sum(Eb_N0_i) / countCycles(Eb_N0_i) / length(EncodeOutput)/q) ...
            ' errF:' num2str(err_Frames(Eb_N0_i))    '   Counts:'  num2str(countCycles(Eb_N0_i))]);
        if err_Frames(Eb_N0_i) >=50 && countCycles(Eb_N0_i)  >=250
            break;
        end
    end
    BER(Eb_N0_i) = err_bit_sum(Eb_N0_i) / (countCycles(Eb_N0_i) * length(EncodeOutput)*q);
    FER(Eb_N0_i) = err_Frames(Eb_N0_i) / countCycles(Eb_N0_i);
    save 6_4.mat
end
save(mfilename,'Eb_N0_dB','BER','FER','countCycles','err_Frames','err_bit_sum');