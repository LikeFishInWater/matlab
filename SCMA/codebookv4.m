%% my codebook
clear;
S=zeros(4,6,4);
C=zeros(4,3);%% constellation
C(1,1)=1+1i;
C(2,1)=-1+1i;
C(3,1)=-1-1i;
C(4,1)=1-1i;
C(1,2)=1/2+1/2*i;
C(2,2)=-1/2+1/2*i;
C(3,2)=-1/2-1/2*i;
C(4,2)=1/2-1/2*i;
C(1,3)=1/4+1/4*i;
C(2,3)=-1/4+1/4*i;
C(3,3)=-1/4-1/4*i;
C(4,3)=1/4-1/4*i;

F=[1 1 0 0;
    1 0 1 0;
    1 0 0 1;
    0 1 1 0;
    0 1 0 1;
    0 0 1 1;]';
h=zeros(4,3);
for m=1:4        
    h(m,:)=find(F(m,:));
end
for ii=1:4
    for j=1:3
        S(ii,h(ii,j),:)=C(:,j);
    end
end
% alpha=[1 2 3 4 5 6;5 2 3 1 4 6;1 4 3 6 5 2;1 2 3 4 5 6;];
% for i=1:4
%     S(i,:,:)=S(i,alpha(i,:),:);
% end
S=S/sqrt(2.6);
save Sv4 S F;

