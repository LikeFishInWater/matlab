clear;

% S=zeros(4,4);
% for m=1:4
%     for n=1:4
% 		S(m,n)=exp(i*pi/4*n);     
%     end
% end
% H=1/sqrt(3)*S;
S=[1 1 ;
    1i 1i ;
    1/3-1/3*i 1/3-1/3*i]';
H=1/sqrt(2.2276)*S;
for ii=1:8
	a=de2bi(ii-1,3);
	a=a*2-1;
    atemp=H*a';
	yr(ii)=atemp(1);
end

x_decode=zeros(3,1);
decode=zeros(3,1);


tic;
EsNo=5:15;
for ir=1:length(EsNo)
    snr=EsNo(ir);
    flag=1;
    fe(ir)=0;
    be(ir)=0;
    total(ir)=0;
    while flag
        total(ir)=total(ir)+1;
        x=2*randi([0 1],3,1)-1;  
        temp=H*x; 
        y=awgn(temp,snr);
        sigma=sqrt(1/2*1/(10^(snr/10)));
        for ii=1:8
			p(ii)=(y(1)-yr(ii))^2+(y(2)-yr(ii))^2;
		end
		[m n]=min(p);
		decode=2*de2bi(n-1,3)-1;

        if decode==x'
            fe(ir)=fe(ir);
        else
            fe(ir)=fe(ir)+1 ;          
        end 
        if(fe(ir)==30)
            flag=0;
        end
        for ii=1:3
            if decode(ii)==x(ii)
                be(ir)=be(ir);
            else
                be(ir)=be(ir)+1;
            end
        end
    end
    ratio_fe(ir)=fe(ir)/total(ir)
    ratio_be(ir)=be(ir)/total(ir)/3
toc;
end
semilogy(EsNo,ratio_be,'bo-');
legend('lds1216');hold on;
xlabel('snr');ylabel('ber');title('lds-cdma12*16');

grid on;



