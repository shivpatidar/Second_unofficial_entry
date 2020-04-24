function [ segs ] = segment_ecg_RtoR( ecg,QRS,fs )

%-------------------FULL SEGMENTATION----------------
tm=QRS;
segs{1,1}=ecg(1:tm(1));
for i=2:length(tm)-2
   
   
    segs{1,i}=ecg(tm(i):tm(i+1)); 

end
segs{1,i+1}=ecg(tm(end):end);
%%plot and see
% for i=1:size(segs,2)
%     plot(segs{1,i})
%     pause(1)
% end

%connecting all segments
% aa1=[];
% for i=1:size(segs,2)
%     aa2=segs{1,i};
%     aa1=[aa1; aa2];
%     
% end


