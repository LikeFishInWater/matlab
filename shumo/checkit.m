function [ sample,training,checking ] = checkit( index,A,Atable )
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
    [m n]=size(index);
    sample=A(:,index');
%     for i=1:1000
%         for j=1:m
%             if sample(i,j)==Atable(index(j,1),index(j,2))
%                 sample(i,j)=posibility(index(j,1),index(j,2));
%             else
%                 sample(i,j)=0;
%             end
%         end
%     end
    for i=1:1000
        for j=1:n
            if sample(i,j)==Atable(1,index(j))
                sample(i,j)=1;
            elseif sample(i,j)==Atable(2,index(j))
                sample(i,j)=2;                    
            else
                sample(i,j)=3;
            end
        end
    end
    training=sample(251:750,:)';
    clear checking;
    checking(1:250,:)=sample(1:250,:);
    checking(251:500,:)=sample(751:1000,:);
    checking=checking';
end

