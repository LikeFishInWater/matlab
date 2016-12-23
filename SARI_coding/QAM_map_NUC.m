function y = QAM_map_NUC(x,constellation)
[N, M]=size(constellation);
msgMD=reshape(x,log2(M),length(x)/log2(M));
msgMD=msgMD.';
[N_x, M_x]=size(msgMD);
for i=1:N_x
    msg(:,i)=bi2de(msgMD(i,:),'left-msb');
end
y=constellation(:,msg+1);