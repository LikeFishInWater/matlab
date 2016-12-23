function LLR=QAM_demap(y,N0,constellation)
%y:信道输出
%N0:噪声方差
[N_y, M_y]=size(y);
[N, M]=size(constellation);
for idx = 1 : log2(M)
    constellation0idx(idx, :) = floor((0 : M / 2 - 1) / 2 ^ (idx - 1)) * 2 ^ (idx - 1) + (1 : M / 2);
end
for idx = 1 : log2(M)
    constellation1idx(idx, :) = ceil((1 : M / 2 ) / 2 ^ (idx - 1)) * 2 ^ (idx - 1) + (1 : M / 2);
end
for idx = 1:log2(M)
    for i=1:N_y
        a = -abs((y(i,:)-constellation(:,constellation0idx(((log2(M)+1)-idx), (1:M/2)))).^2) /N0(i,:);
        b = -abs((y(i,:)-constellation(:,constellation1idx(((log2(M)+1)-idx), (1:M/2)))).^2)/N0(i,:);
        if a>700
            a = 700;
        end
        if a<-700
            a = -700;
        end
        if b>700
            b = 700;
        end
        if b<-700
            b = -700;
        end
        c = sum(exp(a))/sum(exp(b));
        if c>10^200
            c = 10^200;
        end
        if c<1/10^200
            c = 1/10^200;
        end
        LLR(idx,i)=log(c);
    end
end
LLR=reshape(LLR,1,[]).';