clear;clc;
len=1000;

pn=rand(1,len);
sys=tf(1,[2 0 1]);
x=lsim(sys,pn,1:len);
plot(pn);
pause
plot(x)
