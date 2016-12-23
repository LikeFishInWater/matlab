multi_phenos=load('multi_phenos.txt');
multi_sum=sum(multi_phenos,2);
multi_chiP=zeros(1,9445);
r=20;
for j=1:9445
    multi_table=zeros(11+1,3+1);
    for i=1:1000
        multi_table(multi_sum(i)+1,Asimp(i,j))=multi_table(multi_sum(i)+1,Asimp(i,j))+1;
    end
    multi_table(12,:)=sum(multi_table);
    multi_table(:,4)=sum(multi_table,2);
    sumt=0;
    for m=1:11
        for n=1:3
            expect=multi_table(m,4)*multi_table(12,n)/multi_table(12,4);
            sumt=sumt+(multi_table(m,n)-expect)^2/expect;
        end
    end
    multi_chiP(j)=chi2pdf(sumt,r);
end
[I4,B4]=sort(multi_chiP);