function [F1,R1,F2,R2]=makerun(y,M,r);
%function [F1,R1,F2,R2]=makerun(y,M,r);
%
%Input
%
%y input data
%M maximum template length
%r matching tolerance
%
%Output
%
%F1 matches with future points
%R1 runs with future points
%F2 matches with past points
%R2 runs with past points
n=length(y);
run1=zeros(1,n);
MM=2*M;
R1=zeros(n,MM);
R2=zeros(n,MM);
F=zeros(n,M);
F1=zeros(n,M);
for i=1:(n-1)
j=(i+1):n;
match=abs(y(j)-y(i))<r;
k=find(match);
nj=length(j);
run=zeros(1,length(j));
run(k)=run1(k)+1;
for m=1:M
k=find(run>=m);
nm=length(k);
F1(i,m)=nm;
F(i,m)=F(i,m)+nm;
F(i+k,m)=F(i+k,m)+1;
end
nj=min(nj,MM);
k=(1:nj);
R1(i,k)=run(k);
run1=run;
end
for i=1:n
nj=min(MM,i-1);
for j=1:nj
R2(i,j)=R1(i-j,j);
end
end
F2=F-F1;

