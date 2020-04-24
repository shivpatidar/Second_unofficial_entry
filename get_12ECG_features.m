function [fb_v1, fb_avr, avb_feat, lbbb_feat, normal_feat, fb_ar] = get_12ECG_features(data, header_data)

       % addfunction path needed
        addpath(genpath('Tools/'))
        %load('HRVparams_12ECG','HRVparams')
        load('dist_final_temp','template')
        load('data_sig','data_sig')
	% read number of leads, sample frequency and gain from the header.	

	[recording,Total_time,num_leads,Fs,gain,age,sex]=extract_data_from_header(header_data);

    fb_v1=get_fbfeat(data,7);
    fb_avr=get_fbfeat(data,4);
    fb_lead2=get_fbfeat(data,2);
    pr_feat=get_prfeat(data,1);
    dist_feat=get_distfeat(data,7,template);
    ar_feat=get_arfeat(data,8,data_sig);
   avb_feat=[pr_feat fb_v1];
   lbbb_feat=[fb_v1 dist_feat];
   normal_feat=[fb_avr fb_v1 fb_lead2 age];
   fb_ar=[fb_v1 ar_feat];
    
  
end

function [feat] = get_fbfeat(data,lead)
fs=500;
        if size(data,2)<5000
            ecg11 = [data(lead,:) data(lead,:)];
            ecg22=ecg11(:,1:5000);
        else
        ecg22=data(lead,1:5000);    
        end
        ecg=BP_filter_ECG(ecg22,fs);
        feat=feat_29_2020(ecg);
end

function [feat] = get_prfeat(data,lead)
fs=500;
samp_len=250;
PR_seg=[];
piks=[];
PR_idx=[];
P_idx=[];
PR_interval=[];

ecg=data(lead,:);
ecg=BP_filter_ECG(ecg,500);
[QRS,sign,en_thres] = qrs_detect2(ecg',0.25,0.6,fs);%Detecting QRS ( Note: Included as it is from the sample file)
for i=1:1:size(QRS,2)-3
    PR_seg=ecg(1,(QRS(1,i+1)-(0.25*fs)):QRS(1,i+1));
    piks=findpeaks(PR_seg);
    m=max(piks);
    if isempty(m)
      PR_interval(i,1)=60;  
    else
    PR_idx=find(PR_seg==m);
    P_idx=(QRS(1,i+1)-(0.25*fs))+PR_idx(1,1);
    PR_interval(i,1)=QRS(1,i+1)-P_idx;
    end
end
    pr=rmoutliers(PR_interval,'percentile',[20 100]);
    if isempty(pr)
        feat=[60 4 20 60 1 8];
    else
        
    feat=[mean(pr) var(pr) std(pr) median(pr) skewness(pr) kurtosis(pr)];
   
    end
end

function [feat] = get_distfeat(data,lead,feature2)
fs=500;
ar_order=8;
samp_len=250;
ecg=data(lead,:);
[QRS,sign,en_thres] = qrs_detect2(ecg',0.25,0.6,fs);%Detecting QRS ( Note: Included as it is from the sample file)
for i=1:1:size(QRS,2)-3
ecg_seg2(i,:)=ecg((QRS(1,i+1)-(0.2*fs)):(QRS(1,i+1)-(0.2*fs)+samp_len-1));
end
try
for i=1:1:size(ecg_seg2,1)
feature1(i,:) = getarfeat(ecg_seg2(i,:)',ar_order,samp_len,samp_len);
end
k=1;
for i=1:1:20
for j=1:1:size(feature1,1)
pf2=abs(fft((feature1(j,:))).^2);
pf1=abs(fft((feature2(i,:))).^2);
d_itar(k,:) =distitar(feature2(i,:),feature1(j,:),'d');
d_itpf(k,:)=distitpf(pf1,pf2,'d');
% d_eu(k,:)=disteusq(x,y,mode,w);
d_itsar(k,:)=distisar(feature2(i,:),feature1(j,:),'d');
d_copf(k,:)=distchpf(pf1,pf2,'d');
d_coar(k,:)=distchar(feature2(i,:),feature1(j,:),'d');
d_itspf(k,:)=distispf(pf1,pf2,'d');
k=k+1;
end
end
feat=[mean(d_coar) mean(d_itar) mean(d_itsar) mean(d_copf) mean(d_itpf) mean(d_itspf)];
catch
    feat=[1 1 1 1 1 1];
end
end

function [feat] = get_arfeat(data,lead,data_sig)
samp_len=250;
fs=500;
ecg=data(lead,:);
     ecg=BP_filter_ECG(ecg,fs);
[QRS,sign,en_thres] = qrs_detect2(ecg',0.25,0.6,fs);%Detecting QRS ( Note: Included as it is from the sample file)
if size(QRS,2)==3
    for i=1:1:size(QRS,2)-2
        ecg_seg2(i,:)=ecg((QRS(1,i+1)-(0.2*fs)):(QRS(1,i+1)-(0.2*fs)+samp_len-1));
    end
else
for i=1:1:size(QRS,2)-3
        ecg_seg2(i,:)=ecg((QRS(1,i+1)-(0.2*fs)):(QRS(1,i+1)-(0.2*fs)+samp_len-1));
end
end
[m,n]=size(ecg_seg2);
    tdata=reshape(ecg_seg2,1,m*n);
    if m<10
        tdata=[tdata tdata tdata tdata tdata];
        tdata=[tdata tdata];
        end
    if isempty(tdata)
           data=data_sig; 
           [m,n]=size(data);
    tdata=reshape(data,1,m*n);
        end
    
        train=tdata(1,1:2500); 
        
timeWindow = 250;
ARorder = 4;
MODWPTlevel = 4;
[feat,featureindices] = ...
    helperExtractFeatures(train,timeWindow,ARorder,MODWPTlevel);
end

function [filt_signal1] = BP_filter_ECG(ecg,fs)

ecg=ecg;
fs=fs;

d = designfilt('bandpassiir','FilterOrder',6, ...
    'HalfPowerFrequency1',1,'HalfPowerFrequency2',35, ...
    'SampleRate',fs);
%% Filtering 
    filt_signal1=filtfilt(d,ecg);
   
end

function feat = getarfeat(x,order,winsize,wininc,datawin,dispstatus)

if nargin < 6
    if nargin < 5
        if nargin < 4
            if nargin < 3
                winsize = size(x,1);
            end
            wininc = winsize;
        end
        datawin = ones(winsize,1);
    end
    dispstatus = 0;
end

datasize = size(x,1);
%Nsignals = size(x,2);
Nsignals = 1;
numwin = floor((datasize - winsize)/wininc)+1;

% allocate memory
%feat = zeros(numwin,Nsignals*order);
feat = zeros(numwin,order);

if dispstatus
    h = waitbar(0,'Computing AR features...');
end

st = 1;
en = winsize;

for i = 1:numwin
   if dispstatus
       waitbar(i/numwin);
   end
   curwin = x(st:en,:).*repmat(datawin,1,Nsignals);

   cur_xlpc = real(lpc(curwin,order)');
   cur_xlpc = cur_xlpc(2:(order+1),:);
   feat(i,:) = reshape(cur_xlpc,Nsignals*order,1)';
   
   st = st + wininc;
   en = en + wininc;
end

if dispstatus
    close(h)
end
end


function [trainFeatures, featureindices] = helperExtractFeatures(trainData,T,AR_order,level)
% This function is only in support of XpwWaveletMLExample. It may change or
% be removed in a future release.
trainFeatures = [];
%testFeatures = [];

for idx =1:size(trainData,1)
    x = trainData(idx,:);
    x = detrend(x,0);
    arcoefs = blockAR(x,AR_order,T);
    se = shannonEntropy(x,T,level);
    [cp,rh] = leaders(x,T);
    wvar = modwtvar(modwt(x,'db2'),'db2');
    trainFeatures = [trainFeatures; arcoefs se cp rh wvar']; %#ok<AGROW>

end

% for idx =1:size(testData,1)
%     x1 = testData(idx,:);
%     x1 = detrend(x1,0);
%     arcoefs = blockAR(x1,AR_order,T);
%     se = shannonEntropy(x1,T,level);
%     [cp,rh] = leaders(x1,T);
%     wvar = modwtvar(modwt(x1,'db2'),'db2');
%     testFeatures = [testFeatures;arcoefs se cp rh wvar']; %#ok<AGROW>
% 
% end

featureindices = struct();
% 4*8
featureindices.ARfeatures = 1:32;
startidx = 33;
endidx = 33+(16*8)-1;
featureindices.SEfeatures = startidx:endidx;
startidx = endidx+1;
endidx = startidx+7;
featureindices.CP2features = startidx:endidx;
startidx = endidx+1;
endidx = startidx+7;
featureindices.HRfeatures = startidx:endidx;
startidx = endidx+1;
endidx = startidx+13;
featureindices.WVARfeatures = startidx:endidx;
end


function se = shannonEntropy(x,numbuffer,level)
numwindows = numel(x)/numbuffer;
y = buffer(x,numbuffer);
se = zeros(2^level,size(y,2));
for kk = 1:size(y,2)
    wpt = modwpt(y(:,kk),level);
    % Sum across time
    E = sum(wpt.^2,2);
    Pij = wpt.^2./E;
    % The following is eps(1)
    se(:,kk) = -sum(Pij.*log(Pij+eps),2);
end
se = reshape(se,2^level*numwindows,1);
se = se';
end


function arcfs = blockAR(x,order,numbuffer)
numwindows = numel(x)/numbuffer;
y = buffer(x,numbuffer);
arcfs = zeros(order,size(y,2));
for kk = 1:size(y,2)
    artmp =  arburg(y(:,kk),order);
    arcfs(:,kk) = artmp(2:end);
end
arcfs = reshape(arcfs,order*numwindows,1);
arcfs = arcfs';
end


function [cp,rh] = leaders(x,numbuffer)
y = buffer(x,numbuffer);
cp = zeros(1,size(y,2));
rh = zeros(1,size(y,2));
for kk = 1:size(y,2)
    [~,h,cptmp] = dwtleader(y(:,kk));
    cp(kk) = cptmp(2);
    rh(kk) = range(h);
end
end
