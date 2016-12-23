clear;
%indicator matrix
maxIte=1;
%% read F matrix
BinaryFile = 'MatszhN24_12_1_2.txt';
fid=fopen(BinaryFile);
line=str2num(fgetl(fid));
K=line(1);
N=line(2); 
line=str2num(fgetl(fid));
dv=line(1);
dc=line(2);
fgetl(fid);fgetl(fid);
H_row_id = [];
H_col_id = [];
F=zeros(N,K);
for col = 1:K
    row_id = str2num(fgetl(fid));
    F(row_id,col)=1;
end
fclose(fid);
%% 
%spreading signature matrices
h=zeros(N,dc);
for m=1:N         
    h(m,:)=find(F(m,:));
end
v=zeros(dv,K);
for m=1:K
    v(:,m)=find(F(:,m));
end
status=zeros(dc,2^dc);
for t=1:2^dc
    a=de2bi(t-1,dc);
    a=a*2-1;
    status(:,t)=a;
end
S=zeros(N,K);
for m=1:N
    for n=1:K
        if F(m,n)==1
            S(m,n)=sqrt(1/dc)*exp(i*pi/K*n);
        end
        % S(m,h(m,n))=sqrt(1/4)*exp(i*pi/2*n);
    end
end
H=S;

x_decode=zeros(K,1);
decode=zeros(K,1);

tic;
EsNo=0:1:10;
for ir=1:length(EsNo)
    snr=EsNo(ir);
    flag=1;
    fe(ir)=0;
    be(ir)=0;
    total(ir)=0;
    while flag
        total(ir)=total(ir)+1;
        x=2*randi([0 1],K,1)-1;  
        % x=[0 1 0 1 0 0 1 0 1 1 1 0 1 0 1 0]';
        temp=H*x; 
        
        y=awgn(temp,snr);
        sigma=sqrt(1/2*1/(10^(snr/10)));
        %%%%%%%%%%%%%%% MPA decoder %%%%%%%%%%%%%%%%%%%%%
        D=zeros(N,K);
        for j=1:maxIte
            %%%%%%%%%%%%%%%%%Uk to Cn
            for m=1:K;
                he=0;
                for n=1:dv
                    he=he+D(v(n,m),m);
                end
                for n=1:dv
                    D(v(n,m),m)=he-D(v(n,m),m);
                end
            end  
            %%%%%%%%%%%%%Cn to Uk                 
            D_temp=D;
            for n=1:N  
                for c=1:dc
                    rpc=1;rnc=1;
                    for t=1:2^dc
                        if(status(c,t)==1)
                            sum1=0;sum2=0;
                            for tt=1:dc
                                sum1=sum1+status(tt,t)/2*D(n,h(n,tt));
                                sum2=sum2+status(tt,t)*H(n,h(n,tt));
                            end
                            rp(rpc)=sum1-status(c,t)/2*D(n,h(n,c))-1/(2*sigma^2)*(abs(y(n)-sum2))^2;
                            rpc=rpc+1;
                        else
                            sum1=-status(c,t)/2*D(n,h(n,c));sum2=0;
                            for tt=1:dc
                                sum1=sum1+status(tt,t)/2*D(n,h(n,tt));
                                sum2=sum2+status(tt,t)*H(n,h(n,tt));
                            end
                            rn(rnc)=sum1-1/(2*sigma^2)*(abs(y(n)-sum2))^2;
                            rnc=rnc+1;
                        end
                    end
                    % D_temp(n,h(n,c))=log(sum(exp(rp))/sum(exp(rn)));
                    ptem=rp(1);
                    ntem=rn(1);
                    for t=2:2^(dc-1)
                        ptem=max(ptem,rp(t));
                        ntem=max(ntem,rn(t));
                    end
                    D_temp(n,h(n,c))=ptem-ntem;
                end
            end
            D=D_temp;
            %%%%%%%%%%%%%µü´úÖÐ¼ì²â
            for k=1:K
                he(k)=0;
                for n=1:N
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
        for k=1:K
            he(k)=0;
            for n=1:N
                he(k)=he(k)+D(n,k);
            end
            if real(he(k))>0
                x_decode(k)=1;
            else
                x_decode(k)=-1;
            end
        end

        for p=1:K
            if x_decode(p)==x(p)
                be(ir)=be(ir);
            else
                be(ir)=be(ir)+1;
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
    fer(ir)=fe(ir)/total(ir)
    ber(ir)=be(ir)/total(ir)/K
toc;
end
semilogy((1+ir_snr):(ir_max+ir_snr),ber,'bo-');
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


