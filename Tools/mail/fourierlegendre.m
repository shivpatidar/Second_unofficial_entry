function [ a3 ] = fourierlegendre(RR)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%s=(RR-mean(RR))./max(RR); can be try with mean substraction 

s=RR./max(RR);

MM=length(s); %order of FB expansion
%computation of roots of bessel function Jo(x)

x = -1:2/MM:1; 
   
        alfa(1:MM+1,1:MM+1)=legendre(MM,x, 'norm' );
        
        abc1=s*alfa;
        
        
end

