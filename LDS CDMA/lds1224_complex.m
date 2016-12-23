clc;
clear;
%indicator matrix
F=[ -1  -1   1  -1  -1   1  -1  -1  -1   1  -1  -1  -1  -1   1  -1  -1   1  -1  -1  -1  -1  -1   1 
 -1  -1   1  -1   1  -1  -1  -1  -1  -1  -1   1  -1  -1  -1   1  -1  -1   1  -1   1  -1  -1  -1 
 -1  -1  -1   1  -1  -1   1  -1  -1  -1   1   1  -1  -1  -1  -1  -1  -1   1  -1  -1  -1  -1   1 
 -1  -1  -1   1   1  -1  -1  -1   1  -1  -1  -1  -1  -1   1  -1   1  -1  -1  -1   1  -1  -1  -1 
  1  -1  -1  -1  -1  -1  -1   1  -1  -1  -1  -1   1  -1   1  -1  -1  -1   1  -1  -1   1  -1  -1 
  1  -1  -1  -1   1  -1  -1  -1  -1  -1   1  -1  -1  -1  -1   1  -1   1  -1  -1  -1  -1   1  -1 
 -1   1  -1  -1  -1  -1  -1   1  -1  -1  -1   1  -1   1  -1  -1   1  -1  -1  -1  -1  -1   1  -1 
 -1   1  -1  -1  -1   1  -1   1  -1  -1  -1  -1  -1  -1  -1   1  -1  -1  -1   1  -1  -1  -1   1 
  1  -1  -1  -1  -1  -1  -1  -1  -1   1   1  -1  -1   1  -1  -1  -1  -1  -1   1   1  -1  -1  -1 
 -1  -1  -1   1  -1   1  -1  -1   1  -1  -1  -1  -1   1  -1  -1  -1   1  -1  -1  -1   1  -1  -1 
 -1  -1   1  -1  -1  -1   1  -1   1  -1  -1  -1   1  -1  -1  -1  -1  -1  -1   1  -1  -1   1  -1 
 -1   1  -1  -1  -1  -1   1  -1  -1   1  -1  -1   1  -1  -1  -1   1  -1  -1  -1  -1   1  -1  -1 ];
F=(F+1)/2;
% S=zeros(12,24);
% for m=1:24
%     t=1;
%     st=rand(1,3);
%     for n=1:12
%          if F(n,m)==1
%             S(n,m)=sqrt(1/6)*exp(i*2*pi*st(t));
%             t=t+1;
%         end
%     end
% end
S=zeros(12,24);
for m=1:24
    for n=1:12
        if F(n,m)==1
            S(n,m)=sqrt(1/4)*exp(1i*pi/24*m);
        end
    end
end
S=S*sqrt(1/1.5);
%effective received signature matrix
H=S;

x_decode=zeros(24,1);
decode=zeros(24,1);

h=zeros(12,6);
for m=1:12          
    h(m,:)=find(F(m,:));
end
tic;
for ir=1:2
    snr=ir+10;
    flag=1;
    fe(ir)=0;
    be(ir)=0;
    total(ir)=0;
	user_error=zeros(1,24);
    while flag
        total(ir)=total(ir)+1;
        x=2*randi([0 1],24,1)-1;  
        temp=H*x; 
        y=awgn(temp,snr);  
        sigma=sqrt(1/2*1/(10^(snr/10)));
        %%%%%%%%%%%%%%%MPA decoder%%%%%%%%%%%%%%%%%%%%%
        D=zeros(12,24);
        for j=1:15
            %%%%%%%%%%%%%%%%%Uk to Cn
            for m=1:24;
                he=0;
                for n=1:12
                    he=he+D(n,m);
                end
                for n=1:12
                    if(F(n,m))
                        D(n,m)=(he-D(n,m));
                    end
                end
            end        
            %%%%%%%%%%%%%Cn to Uk                 
            D_temp=D;
            for n=1:12  
                for c=1:6
                    rpc=1;rnc=1;
                    for t=1:64
                        a=de2bi(t-1,6);
                        a=a*2-1;
                        if(a(c)==1)
                            rp(rpc)=a(1)/2*D(n,h(n,1))+a(2)/2*D(n,h(n,2))+a(3)/2*D(n,h(n,3))+a(4)/2*D(n,h(n,4))+a(5)/2*D(n,h(n,5))+a(6)/2*D(n,h(n,6))-a(c)/2*D(n,h(n,c))-1/(2*sigma^2)*(abs(y(n)-a(1)*H(n,h(n,1))-a(2)*H(n,h(n,2))-a(3)*H(n,h(n,3))-a(4)*H(n,h(n,4))-a(5)*H(n,h(n,5))-a(6)*H(n,h(n,6))))^2;
                            rpc=rpc+1;
                        else
                            rn(rnc)=a(1)/2*D(n,h(n,1))+a(2)/2*D(n,h(n,2))+a(3)/2*D(n,h(n,3))+a(4)/2*D(n,h(n,4))+a(5)/2*D(n,h(n,5))+a(6)/2*D(n,h(n,6))-a(c)/2*D(n,h(n,c))-1/(2*sigma^2)*(abs(y(n)-a(1)*H(n,h(n,1))-a(2)*H(n,h(n,2))-a(3)*H(n,h(n,3))-a(4)*H(n,h(n,4))-a(5)*H(n,h(n,5))-a(6)*H(n,h(n,6))))^2;
                            rnc=rnc+1;
                        end
                    end
                    D_temp(n,h(n,c))=log(sum(exp(rp))/sum(exp(rn)));
                end
            end
            D=D_temp;
            %%%%%%%%%%%%%µü´úÖÐ¼ì²â
            for k=1:24
                he=0;
                for n=1:12
                    he=he+D(n,k);
                end
                if he>0
                    decode(k)=1;
                else
                    decode(k)=-1;
                end           
            end
            if decode==x
                break;
            end
        end
        for k=1:24
            he=0;
            for n=1:12
                he=he+D(n,k);
            end
            if he>0
                x_decode(k)=1;
            else
                x_decode(k)=-1;
            end
        end

        for p=1:24
            if x_decode(p)==x(p)
                be(ir)=be(ir);
            else
                be(ir)=be(ir)+1;
				user_error(p)=user_error(p);
            end
        end
        if x_decode==x
            fe(ir)=fe(ir);
        else
            fe(ir)=fe(ir)+1;
        end 
        if(fe(ir)==30)
            flag=0;
        end
    end
	worst_user(ir)=max(user_error)/total(ir);
	best_user(ir)=min(user_error)/total(ir);
    ratio_fe(ir)=fe(ir)/total(ir);
    ratio_be(ir)=be(ir)/total(ir)/24
toc;
end
semilogy(1:10,ratio_be,'bo-');
legend('lds1224');
axis([1 10 1e-5 1e0])
grid on;
% str=['figure',num2str(sav)];
% str2=['S_mat',num2str(sav),'.mat'];
% print(1,'-dpng',str);
% save(str2,'S');
save



