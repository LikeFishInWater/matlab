j=1;
for i=1:64
    flag=1;
    t=1;
    while flag && t<100
        if index(i,2)==I(t)
            answer(j)=index(i,2);
            j=j+1;
            flag=0;
        end
        t=t+1;
    end 
end