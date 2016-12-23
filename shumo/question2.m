clear;
load initial.mat;
algorithm='algorithm1';
%% ratio_A health_A total_A
ratio_A=zeros(3,9445);
health_A=zeros(3,9445);
total_A=zeros(3,9445);
for j=1:9445
    for i=1:3
        num_0=0;
        num_1=0;
        for t=1:500
            if A(t,j)==Atable(i,j)
                num_0=num_0+1;
            end
        end
        for t=501:1000
            if A(t,j)==Atable(i,j)
                num_1=num_1+1;
            end
        end
        ratio_A(i,j)=num_0/num_1;
        health_A(i,j)=num_0;
        total_A(i,j)=num_0+num_1;
    end
end
if algorithm=='algorithm1'
    posibility=zeros(3,9445);
    for i=1:3
        for j=1:9445
            posibility(i,j)=posib(total_A(i,j),health_A(i,j));
        end
    end
    t=1;
    clear index;
    for i=1:3
        for j=1:9445
            if posibility(i,j)<0.01
                index(t,:)=[i j posibility(i,j)];
                t=t+1;
            end
        end
    end
    [B,I]=sort(index);
    index0=index(I(:,3),:);
    [m n]=size(index0);
    index1=cell(m,n);
    for i=1:m
        index1{i,2}=AtableRef{index0(i,1),index0(i,2)};
    end
    for i=1:m
        index1{i,1}=Ahead{index0(i,2)};
        index1{i,3}=index0(i,3);
    end
    
elseif algorithm=='algorithm2' %% posibility algorithm2 chi2pdf 
    for j=1:9445
        chi_p(j)=chi_sqr(total_A(:,j),health_A(:,j));
    end
    [B,I]=sort(chi_p);
    % [sample training checking]=checkit(I(1:10),A,Atable);
end
%% common weidian
% j=1;
% for i=1:64
%     flag=1;
%     t=1;
%     while flag && t<100
%         if index(i,2)==I(t)
%             answer(j)=index(i,2);
%             j=j+1;
%             flag=0;
%         end
%         t=t+1;
%     end 
% end
%% sample checking training
% [m n]=size(index);
% sample=A(:,index(:,2)');
% for i=1:1000
%     for j=1:m
%         if sample(i,j)==Atable(index(j,1),index(j,2))
%             sample(i,j)=posibility(index(j,1),index(j,2));
%         else
%             sample(i,j)=0.05;
%         end
%     end
% end
% training=sample(251:750,:)';
% clear checking;
% checking(1:250,:)=sample(1:250,:);
% checking(251:500,:)=sample(751:1000,:);
% checking=checking';