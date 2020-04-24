function [ a3 ] = fourierbessel(RR)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

s=RR;

MM=length(s); %order of FB expansion
%computation of roots of bessel function Jo(x)
if exist('alfa') == 0
    x=2;
    alfa=zeros(1,MM);
    for i=1:MM
        ex=1;
        while abs(ex)>.00001
            ex=-besselj(0,x)/besselj(1,x);
            x=x-ex;
        end
        alfa(i)=x;

        x=x+pi;
    end
end


N=length(s);
nb=1:N;


a=N;

for m1=1:MM
    a3(m1)=(2/(a^2*(besselj(1,alfa(m1))).^2))*sum(nb.*s.*besselj(0,alfa(m1)/a*nb));
end



end

