function [ chi_p ] = chi_sqr(  total_A,health_A)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
    r=2;
    t=0;
    for i=1:3
        t=t+(health_A(i)-total_A(i)/2)^2/(total_A(i)/2)+(total_A(i)-health_A(i)-total_A(i)/2)^2/(total_A(i)/2);
    end
    chi_p=chi2pdf(t,r);
end

