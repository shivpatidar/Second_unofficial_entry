function [ data ] = morphoroughness(segs )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


for i=2:(size(segs,2)-1)
    Kurto=kurtosis(segs{1,i}+10);
    Medn=median(segs{1,i}+10);
    idata(i,1:2)=[  Kurto  Medn ];

end

data1 =  mean(idata(:,1),1);
data2 = median(idata(:,:),1);
data3 =  skewness(idata(:,1),1);

data=[data1 data2 data3 ];
end