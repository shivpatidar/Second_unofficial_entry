function [ data ] = peaksinfo_ecg( ecg,QRS,fs )

tm=QRS;
segs{1,1}=ecg(1:tm(1));
for i=2:length(tm)-2,
   
   
    segs{1,i}=ecg(tm(i):tm(i+1)); 

end
segs{1,i+1}=ecg(tm(end):end);

%plot and see
%figure
idata=[];
for i=2:size(segs,2)-1
%     plot(segs{1,i})
[pks,locs,w,p]=findpeaks(abs(segs{1,i}./max(segs{1,i})),'MinPeakDistance',70,'NPeaks',2,'MinPeakProminence',0.08);
%pause(2)
% if length(pks)<=1
%     [pks,locs,w,p]=findpeaks(segs{1,i}./max(segs{1,i}),'MinPeakDistance',10,'NPeaks',2,'MinPeakProminence',0.05);
% else
% 
% end
% 
% if length(pks)<=1
%     [pks,locs,w,p]=findpeaks(segs{1,i}./max(segs{1,i}),'MinPeakDistance',10,'NPeaks',2,'MinPeakProminence',0.03);
% else
% 
% end
% if length(pks)<=1
%     [pks,locs,w,p]=findpeaks(segs{1,i}./max(segs{1,i}),'MinPeakDistance',10,'NPeaks',2,'MinPeakProminence',0.01);
% else
% 
% end
if length(pks)==1
    pks(2)=0;locs(2)=0;w(2)=0;p(2)=0;
    idata = [idata;pks,locs,w,p];
elseif isempty(pks)==1
    idata=[idata;[0 0 0 0 0 0 0 0]];
else
    
idata = [idata;pks',locs',w',p'];
end

end
data1=[mode(idata,1)];

data=(data1(:,3).*data1(:,1))./( data1(:,4).*data1(:,2)-data1(:,3).*data1(:,1));
end

%data=[mean(sz) median(sz) mode(sz) var(sz) mean(pkspow) median(pkspow) mode(pkspow) var(pkspow) mean(pkspow1) median(pkspow1) mode(pkspow1) var(pkspow1) mean(pkspowwt) median(pkspowwt) mode(pkspowwt) var(pkspowwt) mean(pkspowwt1) median(pkspowwt1) mode(pkspowwt1) var(pkspowwt1) ];


% if isempty(pks)==1 || isempty(pks)==1 
% else
% sz(i)=length(pks);
% pkspow(i)=sum(pks)./max(pks);
% pkspow1(i)=sum(pks./max(pks));
% pkspowwt(i) =sum(pks.*locs)./max(pks.*locs);
% pkspowwt1(i) =sum((pks.*locs)./max(pks.*locs));
%    
% %pause(1)
% end

