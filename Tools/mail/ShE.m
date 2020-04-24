function Hsh = ShE(x)
%shannon entropy
%   Detailed explanation goes here
%x is the input sequence-'coloumn vector'
x1=x'./max(x);
N=length(x1);%number of samples in x
k=N*.1;%k/N=.01 k bins
k1=(max(x1)-min(x1))/k;
y=min(x1):k1:max(x1);%amplitde levels
n=histc(x,y);%Occurence of amplitudes
%plot(y,n)%Histogram plot
for i=1:1:length(y),
Hshp(i)=n(i)*log(n(i));
end
Hshf=-nansum(Hshp);
Hsh=Hshf/log(k);%dont use it for k=1 for it gives inf

end

