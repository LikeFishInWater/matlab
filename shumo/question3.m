%% 
clear;
load initial.mat;
load question2_A2.mat;
scale=100;
index=I(1:scale);
clear ind ind0 ind1 ind2;
ind0=[];
tt=0;
for i=1:300
    [m n]=size(gene{i});
    find=0;
    ind=[];
    for j=1:n
        t=1;
        flag=1;       
        while t<=scale && flag
            if gene{i}(j)==index(t)
                flag=0;
                find=1;
                ind=[ind,j];
            end
            t=t+1;
        end
    end
    if find==1
        ind0=[ind0,i];
        tt=tt+1;
        ind1{tt}=ind;
    end
end
[m n]=size(ind0);
for i=1:n
    ind2{i}=gene{ind0(i)}(ind1{i});
end
for i=1:n
    [mm nn]=size(ind2{i});
    if nn==1
        flag=1;
        t=1;
        while flag
            if ind2{i}==I(t)
                flag=0;
                ind3(i)=B(t);
            end
            t=t+1;
        end
    else 
        ind3(i)=correlationn2(ind2{i},Asimp);
    end
end
[BB II]=sort(ind3);
ind4=ind0(II);
ind5=ind3(II);
ind6=ind2(II);
fid= fopen('gene_bad.txt','wt');
for i=1:n
    tempv=ind2{i};
    [mm,nn]=size(tempv);
    fprintf(fid,'gene_%d: ',ind0(i));
    for j=1:nn       
        fprintf(fid,'%d ',tempv(j));
    end
    % fprintf(fid,'                  %d',ind3(i));
    fprintf(fid,'\n');
end
fclose(fid);
fid=fopen('gene_bad_sort.txt','wt');
for i=1:n
    tempv=ind6{i};
    [mm nn]=size(tempv);
    fprintf(fid,'gene_%d ',ind4(i));
    for j=1:nn
        fprintf(fid,'  %d',ind6{i}(j));
    end
    
    fprintf(fid,'                   %d\n',ind5(i));
end
fclose(fid);