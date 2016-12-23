%%codebook constellation

load Sv3.mat;
for i=1:64
    a(1,i)=mod(i-1,4);
    a(2,i)=mod((i-1-a(1,i))/4,4);
    a(3,i)=((i-1-a(1,i))/4-a(2,i))/4;
end
a=a+1;
Cons=zeros(4,64);
h=zeros(4,3);
for m=1:4        
    h(m,:)=find(F(m,:));
end
for i=1:4
    for j=1:64
        Cons(i,j)=S(i,h(i,1),a(1,j))+S(i,h(i,2),a(2,j))+S(i,h(i,3),a(3,j));
    end
end
figure(1);
scatter(real(Cons(1,:)),imag(Cons(1,:)));
for i=1:64
    text(real(Cons(1,i)),imag(Cons(1,i)),[num2str(a(1,i)),num2str(a(2,i)),num2str(a(3,i))]);
end
figure(2);
scatter(real(Cons(2,:)),imag(Cons(2,:)));
figure(3);
scatter(real(Cons(3,:)),imag(Cons(3,:)));
figure(4);
scatter(real(Cons(4,:)),imag(Cons(4,:)));
for i=1:64
    text(real(Cons(1,i)),imag(Cons(1,i)),[num2str(a(1,i)),num2str(a(2,i)),num2str(a(3,i))]);
end