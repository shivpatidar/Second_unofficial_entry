function Esh = SpE(x)
%spectral entropy
%   Detailed explanation goes here
%x is the input sequence-'coloumn vector'
x1=x'./max(x);
h=spectrum.burg;
xpsd=psd(h,x1);
%fl,fh are min and max frequency limits 0 to 1000Hz here
psd1=xpsd.data./sum(xpsd.data);%normalization such that sum(P)=1
Nf=length(psd1);
for i=1:1:Nf,
Eshp(i)=psd1(i)*log(psd1(i));
end
Eshf=-nansum(Eshp);
Esh=Eshf/log(Nf);

end

