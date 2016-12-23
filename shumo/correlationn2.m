function [ p ] = correlationn2( ind,Asimp)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
    [~,n]=size(ind);
    table=zeros(2,3^n);
    r=3^n-1;
    for i=1:500
        temp=Asimp(i,ind)-1;
        sumt=1;
        for j=1:n
            sumt=sumt+3^(j-1)*temp(j);
        end
        table(1,sumt)=table(1,sumt)+1;       
    end
    for i=501:1000
        temp=Asimp(i,ind)-1;
        sumt=1;
        for j=1:n
            sumt=sumt+3^(j-1)*temp(j);
        end
        table(2,sumt)=table(2,sumt)+1;          
    end
    numt=0;
    for i=1:2
        for j=1:3^n
            if (table(1,j)+table(2,j))~=0
                numt=numt+(table(i,j)-(table(1,j)+table(2,j))/2)^2/((table(1,j)+table(2,j))/2);
            end
        end
    end
    p=chi2pdf(numt,r);
end

