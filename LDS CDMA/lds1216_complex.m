clear;
%indicator matrix
F=[0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0;
    0,0,0,1,0,0,1,0,0,0,0,0,1,0,1,0;
    0,1,0,0,1,0,1,0,0,1,0,0,0,0,0,0;
    0,0,0,1,0,1,0,0,0,0,1,0,0,1,0,0;
    0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,1;
    1,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0;
    0,0,0,1,0,0,0,1,1,0,0,1,0,0,0,0;
    0,1,0,0,0,1,0,0,0,0,0,1,0,0,0,1;
    1,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1;
    0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0;
    0,1,0,0,0,0,0,0,1,0,1,0,1,0,0,0;
    1,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0;];
%spreading signature matrices
h=zeros(12,4);
for m=1:12          
    h(m,:)=find(F(m,:));
end
v=zeros(3,16);
for m=1:16
    v(:,m)=find(F(:,m));
end
status=zeros(4,16);
for t=1:16
    a=de2bi(t-1,4);
    a=a*2-1;
    status(:,t)=a;
end
S=zeros(12,16);
for m=1:12
    for n=1:16
        if F(m,n)==1
            S(m,n)=sqrt(1/4)*exp(i*pi/16*n);
        end
        % S(m,h(m,n))=sqrt(1/4)*exp(i*pi/2*n);
    end
end
H=S;

x_decode=zeros(16,1);
decode=zeros(16,1);


tic;
EsNo=1:1:10;
for ir=1:length(EsNo)
    snr=EsNo(ir);
    flag=1;
    fe(ir)=0;
    be(ir)=0;
    total(ir)=0;
    user_error=zeros(1,16);
    while flag
        total(ir)=total(ir)+1;
        x=2*randi([0 1],16,1)-1;  
        % x=[0 1 0 1 0 0 1 0 1 1 1 0 1 0 1 0]';
        temp=H*x; 
        
        y=awgn(temp,snr);
        sigma=sqrt(1/2*1/(10^(snr/10)));
        %%%%%%%%%%%%%%% MPA decoder %%%%%%%%%%%%%%%%%%%%%
        D=zeros(12,16);
        for j=1:7
            %%%%%%%%%%%%%%%%%Uk to Cn
            for m=1:16;
                he=0;
                for n=1:3
                    he=he+D(v(n,m),m);
                end
                for n=1:3
                    D(v(n,m),m)=he-D(v(n,m),m);
                end
            end  
            %%%%%%%%%%%%%Cn to Uk                 
            D_temp=D;
            for n=1:12  
                for c=1:4
                    rpc=1;rnc=1;
                    for t=1:16
                        if(status(c,t)==1)
                            rp(rpc)=status(1,t)/2*D(n,h(n,1))+status(2,t)/2*D(n,h(n,2))+status(3,t)/2*D(n,h(n,3))+status(4,t)/2*D(n,h(n,4))-status(c,t)/2*D(n,h(n,c))-1/(2*sigma^2)*(abs(y(n)-status(1,t)*H(n,h(n,1))-status(2,t)*H(n,h(n,2))-status(3,t)*H(n,h(n,3))-status(4,t)*H(n,h(n,4))))^2;
                            rpc=rpc+1;
                        else
                            rn(rnc)=status(1,t)/2*D(n,h(n,1))+status(2,t)/2*D(n,h(n,2))+status(3,t)/2*D(n,h(n,3))+status(4,t)/2*D(n,h(n,4))-status(c,t)/2*D(n,h(n,c))-1/(2*sigma^2)*(abs(y(n)-status(1,t)*H(n,h(n,1))-status(2,t)*H(n,h(n,2))-status(3,t)*H(n,h(n,3))-status(4,t)*H(n,h(n,4))))^2;
                            rnc=rnc+1;
                        end
                    end
                    % D_temp(n,h(n,c))=log(sum(exp(rp))/sum(exp(rn)));
                    ptem=rp(1);
                    ntem=rn(1);
                    for t=2:8
                        ptem=max(ptem,rp(t));
                        ntem=max(ntem,rn(t));
                    end
                    D_temp(n,h(n,c))=ptem-ntem;
                end
            end
            D=D_temp;
            %%%%%%%%%%%%%µü´úÖÐ¼ì²â
            for k=1:16
                he(k)=0;
                for n=1:12
                    he(k)=he(k)+D(n,k);
                end
                if real(he(k))>0
                    decode(k)=1;
                else
                    decode(k)=-1;
                end           
            end
            if decode==x
                break;
            end
        end
        for k=1:16
            he(k)=0;
            for n=1:12
                he(k)=he(k)+D(n,k);
            end
            if real(he(k))>0
                x_decode(k)=1;
            else
                x_decode(k)=-1;
            end
        end

        for p=1:16
            if x_decode(p)==x(p)
                be(ir)=be(ir);
            else
                be(ir)=be(ir)+1;
                user_error(p)=user_error(p)+1;
            end
        end
        if x_decode==x
            fe(ir)=fe(ir);
        else
            fe(ir)=fe(ir)+1           
        end 
        if(fe(ir)==30)
            flag=0;
        end
    end
    ratio_fe(ir)=fe(ir)/total(ir)
    ratio_be(ir)=be(ir)/total(ir)/16
toc;
end
semilogy((1+ir_snr):(ir_max+ir_snr),ratio_be,'bo-');
legend('lds1216');hold on;
xlabel('snr');ylabel('ber');title('lds-cdma12*16');
axis([(1+ir_snr) (ir_max+ir_snr) 1e-5 1e-1]);
grid on;
% hold off;
% str=['figure',num2str(sav)];
% str2=['S_mat',num2str(sav),'.mat'];
% print(1,'-dpng',str);
% save(str2,'S');
% end
% ber(sav)=ratio_be;
% end
save lds1216_complex.mat;


