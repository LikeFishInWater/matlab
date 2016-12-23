Atemp=importdata('genotype.dat');%% 源数据有错误
%% code
code={'AA','AT','AC','AG','TA','TT','TC','TG','CA','CT','CC','CG','GA','GT','GC','GG'};
%% Ahead
t=1;
for j=1:9444
    i=t;
    while Atemp{1,1}(1,i)~=' ';
        i=i+1;
    end
    Ahead(j)={Atemp{1,1}(1,t:i-1)};
    t=i+1;
end
Ahead(9445)={'rs7545865'};
%% A
A=zeros(1000,9445);
for i=1:1000
    for j=1:9445
        switch Atemp{i+1,1}(1,(3*(j-1)+1):(3*(j-1)+2)) 
            case 'AA',A(i,j)=1;
            case 'AT',A(i,j)=2;
            case 'AC',A(i,j)=3;
            case 'AG',A(i,j)=4;
            case 'TA',A(i,j)=5;
            case 'TT',A(i,j)=6;
            case 'TC',A(i,j)=7;
            case 'TG',A(i,j)=8;
            case 'CA',A(i,j)=9;
            case 'CT',A(i,j)=10;
            case 'CC',A(i,j)=11;
            case 'CG',A(i,j)=12;
            case 'GA',A(i,j)=13;
            case 'GT',A(i,j)=14;
            case 'GC',A(i,j)=15;
            case 'GG',A(i,j)=16;
        end            
    end
end
for i=1:1000
    switch Atemp{i+1,1}(25417:25418)
        case 'II',A(i,8473)=11;
        case 'ID',A(i,8473)=10;
        case 'DD',A(i,8473)=6;
    end
    switch Atemp{i+1,1}(18208:18209)
        case 'II',A(i,6070)=11;
        case 'ID',A(i,6070)=10;
        case 'DD',A(i,6070)=6;
    end
end
%% Aref
Aref=cell(1000,9445);
for i=1:1000
    for j=1:9445
        Aref{i,j}=code{A(i,j)};
    end
end
%% Atable AtableRef
Atable=zeros(3,9445);
AtableRef=cell(3,9445);
for j=1:9445
    t=1;
    Atable(1,j)=A(t,j);
    AtableRef{1,j}=Aref{t,j};
    t=t+1;
    i=2;
    while i<=3 && t<=1000;
       while ((A(t,j)==Atable(1,j))||(A(t,j)==Atable(2,j)))&&t<1000
           t=t+1;
       end
       Atable(i,j)=A(t,j);
       AtableRef{i,j}=Aref{t,j};
       i=i+1;
    end
end
%% Asimp
Asimp=A;
for i=1:1000
    for j=1:9445
        switch A(i,j)
            case Atable(1,j),Asimp(i,j)=1;
            case Atable(2,j),Asimp(i,j)=2;
            case Atable(3,j),Asimp(i,j)=3;
            otherwise Asimp(i,j)=0;
        end
    end
end
%% flag
clear temp;
for j=1:9445
    for i=1:1000
        temp(i,:)=code{A(i,j)};
    end
    num=0;
    for i=1:1000
        if temp(i,1)==temp(1,1)
            num=num+1;
        end
        if temp(i,2)==temp(1,1)
            num=num+1;
        end 
    end
    if num>100 && num<1900
        flag(j)=1;
    else 
        flag(j)=0;
    end
end
%% gene
for i=1:300
    temp=importdata(['gene_',num2str(i),'.dat']);
    [s t]=size(temp);
    t=1;
    clear numtemp;
    for j=1:s
        flag=1;
        while flag
            if size(temp{j,1})==size(Ahead{1,t})
                if temp{j,1}==Ahead{1,t}
                    numtemp(j)=t;
                    flag=0;
                end
            end
            t=t+1;
            if t>9445
                flag=0;
            end
        end
    end 
    gene(i)={numtemp};
end
%% 

