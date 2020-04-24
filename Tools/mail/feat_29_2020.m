function data = feat_29_2020(ecg)
%
% Sample entry for the 2017 PhysioNet/CinC Challenge.
%
% INPUTS:
% recordName: string specifying the record name to process
%
% OUTPUTS:
% classifyResult: integer value where
%                     N = normal rhythm
%                     A = AF
%                     O = other rhythm
%                     ~ = noisy recording (poor signal quality)
%
% To run your entry on the entire training set in a format that is
% compatible with PhysioNet's scoring enviroment, run the script
% generateValidationSet.m
%
% The challenge function requires that you have downloaded the challenge
% data 'training_set' in a subdirectory of the current directory.
%    http://physionet.org/physiobank/database/challenge/2017/
%
% This dataset is used by the generateValidationSet.m script to create
% the annotations on your training set that will be used to verify that
% your entry works properly in the PhysioNet testing environment.
%
%
% Version 1.0
%
%
% Written by: Shivnarayan Patidar and Ashish Sharma April 2017
%             shivnarayan.patidar@nitgoa.ac.in
%
% Last modified by:
%% rePE.m - Copyright (c) 2016, Valentina Unakafova
%
%load tenthentrymodel1
%load tenthentrymodel2
%load current_model_23_may
%load model_for_29
%load newmodelwithtest
%classifyResult = 'N'; % default output normal rhythm

%% Class determination Normal(N)/AF(A)/Others(O)/Noise(~)

ecg=ecg; %Reading 
fs=500;
%ecg=(ecg-nanmean(ecg))./nanstd(ecg);
[QRS,sign,en_thres] = qrs_detect2(ecg',0.25,0.6,fs);%Detecting QRS ( Note: Included as it is from the sample file)


 if length(QRS)<6
     
     data=ones(1,29);
     
 else
     
%     
%     SHEe=ShE(ecg);%shanon entropy%---------
%     TKurtoecg=kurtosis(Teager(ecg));
%     datafe=[TKurtoecg log(abs(SHEe))];
%     [MD1, re_PE1] = rePE(ecg, 1, 2, 3, 0.01, 0.045);
%     datarepe=[sum(MD1./length(ecg)) ];
%     data=[ datarepe datafe ];
%    
% label = predict(Mdl1,data);%bag of tree based decision
% switch label
%     case 1
%         classifyResult ='N';
%     case 2
%        classifyResult ='A';
%     case 3
%        classifyResult ='O';
%     otherwise
%         classifyResult ='~';
% end
% 
% else
    RR=diff(QRS')/fs;
    if length(RR)<21
        RR=[RR' RR' RR']';
    else
    end
    HR=(1./RR).*60;
    %--------Computing Fourier Bessels of RR and HR interval and respective features--------------
    [ a3hr ] = fourierbessel(HR' );
    [ a3rr ] = fourierbessel(RR' );
    [xfrr,frr] = fft_freq(a3rr,1,[],320);
    % ------------------Direct features-------
    %ecg1=ecg(750:end);
    ecg1=ecg;
    SPEe=SpE(ecg1);%spectral entropy
    SHEe=ShE(ecg1);%shanon entropy%---------
    %--------Teager energy----
    TeagerHR=sum(Teager(HR).^2);
    TKurtoHR=kurtosis(Teager(HR));
    TKurtoecg=kurtosis(Teager(ecg1));
    TSkewa3hr=skewness(Teager(a3hr));
    TStdHR=std(Teager(HR));
   
 
    %-------------Other recent features--------
    SPE=SpE(HR);%Based on HR
   
     SnpHR=sampen(HR,[],0.1,[],0,[]);
     SnpRR=sampen(RR,2,0.1,[],0,[]);

    if (SnpHR(2,1)==Inf)
        SnpHR(2,1)=1;
    else
    end
    
    datafe=[SnpHR(1,1)' SnpRR(2,1)'   ...
        TeagerHR, TKurtoHR,TKurtoecg,TSkewa3hr,TStdHR,...
        SPE, log(abs(SPEe)),log(abs(SHEe))];
    %-------other entropies-------------
    
   
    e8=wentropy(a3rr,'norm',4);
    e11=wentropy(Teager(RR),'sure',0.06);
    e12=wentropy(Teager(RR),'norm',3);
    
    datawen=[ e8   e11 e12 ];
    
    %------rePE-based---features-----------
   
    [MD2, re_PE2] = rePE(RR, 1, 2, 3, 0.01, 0.82);
    [MD4, re_PE4] = rePE(a3hr, 2, 3, 8, 1,141);
    
    datarepe=[ sum(re_PE2.^2) sum(MD2./length(RR))...
         sum(MD4./length(a3hr))];
    
    
    %------------------statistics of RR and HR--------
    
    
    RR_median = median(RR);
    RR_Kurto = kurtosis(RR);
    RR_skew = skewness(RR);
    
    HR_mean = mean(HR);
    HR_median = median(HR);
    HR_mode = mode(HR);
   
    
    datastat = [  RR_median  RR_Kurto RR_skew HR_mean HR_median HR_mode ];
    data1=[datawen datarepe datafe  datastat];
    %---segmentation based morphological features
    
    
    [ segs ] = segment_ecg_RtoR( ecg,QRS,fs );
    data21 = morphoroughness(segs );
    data2=[ data21([2 4])];
    
    datai= [data1 data2 ];
    datam= morpho_ecg2( ecg,QRS );
    
    datam2= morpho_ecg3( ecg,QRS );

%           Snpecg=sampen(ecg(750:end),30,0.1,[],0,[]);
           Snpecg=sampen(ecg,30,0.1,[],0,[]); 

     data5=Snpecg';

      

        dataa=[ datai data5(1,[28 30]) datam(1,[ 6   43  ]) datam2(1,[1 ]) ];
   
    for i2=1:29
        if (isnan(dataa(i2))== 1)
            dataa(i2)=0;
            
        else
        end
        
        if (dataa(i2)==Inf)
            dataa(i2)=1;
            
        else
        end
    end
        
    
    
%zz=[77 -740	807	9 500 285 -397 981 -194	94 -326	648	684	-810 255 -806 -649 -36 902	564	-916 -41 -371 840 -178 -812	550	17	-327 583 315 788 -562 -235	-998 623 -196 403 -434 454	594	801	-10	431	-62	-245 414];
    
data= dataa;
    
 end
   
end


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

%
% FFT_FREQ Compute the FFT and provide the corresponding frequency values.
%
% [xf,f] = fft_freq(x,fs,shift,Nfft)
%

% This function is to get the FFT along with the corresponding frequency
% values. Users have the option of centering the spectrum using (using
% fftshift).
%
% Inputs
%    x:        signal (column vector)
%    fs:       sampling frequency (default = 2*pi)
%    shift:    using fftshift or not (default = false)
%    Nfft:     N-point FFT
%
% Outputs
%    xf:       FFT of x
%    f:        frequency values
%

function [xf,f] = fft_freq(x,fs,shift,Nfft)

N = length(x);

if nargin < 4
    Nfft = N;
    if nargin < 3
        shift = false;
        if nargin < 2
            fs = 2*pi;
        end
    end
end

xf = fft(x, Nfft);

% calculate frequency spacing
df = fs / Nfft;

% calculate unshifted frequency vector
f = (0:(Nfft-1))'*df;

% move all frequencies that are greater than fs/2 to the negative side of the axis
if shift
    xf = fftshift(xf);
    f(f >= fs/2) = f(f >= fs/2) - fs;
    f = fftshift(f);
end
end
function [idata ] = morpho_ecg2( ecg,QRS )

for i=1:(length(QRS)-1)
    
 if QRS(1,i)>20 
     
[aa bb cc dd]= (findpeaks((ecg((QRS(1,i)-20):(QRS(1,i)+50))),'Annotate','extents','SortStr','descend'));


[aa1 bb1 cc1 dd1] = findpeaks((-ecg((QRS(1,i)-20):(QRS(1,i)+50))),'Annotate', 'extents', 'SortStr','descend');

if isempty(cc==1)
    cc=1;
else
end

if  isempty(cc1==1)
    cc1=1;
else
end
if    isempty(dd==1)
    dd=1;
else
end

if   isempty(dd1==1)
    dd1=1;
else
end
   
    
ee = dd(1)/dd1(1);

ff = cc(1)/cc1(1);

gg = (dd(1).*cc(1))/(dd1(1).*cc1(1));

data(i,1:7) = [cc(1) dd(1) cc1(1) dd1(1) ee(1) ff(1) gg(1)];
 else
     
 end
 
end
data1 =  mean(data(:,:),1);
data2 = median(data(:,:),1);
data3 =  mode(data(:,:),1);
data4 = kurtosis(data(:,:),1);
data5 =  skewness(data(:,:),1);
data6 = std(data(:,:),1);
data7 = var(data(:,:),1);


idata=[data1 data2 data3 data4 data5 data6 data7];
    
end
function [data ] = morpho_ecg3( ecg,QRS )

for i=1:(length(QRS)-1)
    
    if QRS(1,i)>20
        
        [aa bb cc dd]= (findpeaks((ecg((QRS(1,i)-20):(QRS(1,i)+50))),'Annotate','extents','SortStr','descend'));
        
        
        [aa1 bb1 cc1 dd1] = findpeaks((-ecg((QRS(1,i)-20):(QRS(1,i)+50))),'Annotate', 'extents', 'SortStr','descend');
        
        if isempty(aa==1)
            aa=1;
        else
        end
        if isempty(aa1==1)
            aa1=1;
        else
        end
        
        if isempty(bb==1)
            bb=1;
        else
        end
        if isempty(bb1==1)
            bb1=1;
        else
        end
        
        if isempty(cc==1)
            cc=1;
        else
        end
        
        if  isempty(cc1==1)
            cc1=1;
        else
        end
        if    isempty(dd==1)
            dd=1;
        else
        end
        
        if   isempty(dd1==1)
            dd1=1;
        else
        end
        
        
        
        
        idata(i,1:6) = [(aa(1)-aa1(1)) (bb(1)-bb1(1)) (cc(1)-cc1(1)) (dd(1)-dd1(1)) aa(1)/aa1(1)...
            bb(1)/bb1(1) ];
    else
        
    end
    
end
data1 =  mean(idata(:,:),1);
data2 = median(idata(:,:),1);
data3 =  mode(idata(:,:),1);
data4 = kurtosis(idata(:,:),1);
data5 =  skewness(idata(:,:),1);
data6 = std(idata(:,:),1);
data7 = var(idata(:,:),1);


data=[data1 data2 data3 data4 data5 data6 data7];
    
end

function [qrs_pos,sign,en_thres] = qrs_detect2(ecg,varargin)
% QRS detector based on the P&T method. This is an offline implementation
% of the detector.
%
% inputs
%   ecg:            one ecg channel on which to run the detector (required)
%                   in [mV]
%   varargin
%       THRES:      energy threshold of the detector (default: 0.6) 
%                   [arbitrary units]
%       REF_PERIOD: refractory period in sec between two R-peaks (default: 0.250)
%                   in [ms]
%       fs:         sampling frequency (default: 1KHz) [Hz]
%       fid_vec:    if some subsegments should not be used for finding the
%                   optimal threshold of the P&Tthen input the indices of
%                   the corresponding points here
%       SIGN_FORCE: force sign of peaks (positive value/negative value).
%                   Particularly usefull if we do window by window detection and want to
%                   unsure the sign of the peaks to be the same accross
%                   windows (which is necessary to build an FECG template)
%       debug:      1: plot to bebug, 0: do not plot
%
% outputs
%   qrs_pos:        indexes of detected peaks (in samples)
%   sign:           sign of the peaks (a pos or neg number)
%   en_thres:       energy threshold used
%
%
%
% Physionet Challenge 2014, version 1.0
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar
% Oxford university, Intelligent Patient Monitoring Group
% joachim.behar@eng.ox.ac.uk
%
% Last updated : 13-09-2014
% - bug on refrac period fixed
% - sombrero hat for prefiltering added
% - code a bit more tidy
% - condition added on flatline detection for overall segment (if flatline 
% then returns empty matrices rather than some random stuff)
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

% == managing inputs
WIN_SAMP_SZ = 7;
REF_PERIOD = 0.250; 
THRES = 0.6; 
fs = 300; 
fid_vec = [];
SIGN_FORCE = [];
debug = 0;

switch nargin
    case 1
        % do nothing
    case 2
        REF_PERIOD=varargin{1};
    case 3
        REF_PERIOD=varargin{1}; 
        THRES=varargin{2};
    case 4
        REF_PERIOD=varargin{1}; 
        THRES=varargin{2};
        fs=varargin{3};  
    case 5
        REF_PERIOD=varargin{1}; 
        THRES=varargin{2}; 
        fs=varargin{3}; 
        fid_vec=varargin{4};
    case 6
        REF_PERIOD=varargin{1}; 
        THRES=varargin{2}; 
        fs=varargin{3};
        fid_vec=varargin{4}; 
        SIGN_FORCE=varargin{5};
    case 7
        REF_PERIOD=varargin{1}; 
        THRES=varargin{2}; 
        fs=varargin{3};
        fid_vec=varargin{4};
        SIGN_FORCE=varargin{5};          
        debug=varargin{6};
    case 8
        REF_PERIOD=varargin{1}; 
        THRES=varargin{2}; 
        fs=varargin{3};
        fid_vec=varargin{4};
        SIGN_FORCE=varargin{5};          
        debug=varargin{6};
        WIN_SAMP_SZ = varargin{7};
    otherwise
        error('qrs_detect: wrong number of input arguments \n');
end


[a b] = size(ecg);
if(a>b); NB_SAMP=a; elseif(b>a); NB_SAMP=b; ecg=ecg'; end;
tm = 1/fs:1/fs:ceil(NB_SAMP/fs);

% == constants
MED_SMOOTH_NB_COEFF = round(fs/100);
INT_NB_COEFF = round(WIN_SAMP_SZ*fs/256); % length is 30 for fs=256Hz  
SEARCH_BACK = 1; % perform search back (FIXME: should be in function param)
MAX_FORCE = []; % if you want to force the energy threshold value (FIXME: should be in function param)
MIN_AMP = 0.1; % if the median of the filtered ECG is inferior to MINAMP then it is likely to be a flatline
               % note the importance of the units here for the ECG (mV) 
NB_SAMP = length(ecg); % number of input samples

try
    % == Bandpass filtering for ECG signal
    % this sombrero hat has shown to give slightly better results than a
    % standard band-pass filter. Plot the frequency response to convince
    % yourself of what it does
    b1 = [-7.757327341237223e-05  -2.357742589814283e-04 -6.689305101192819e-04 -0.001770119249103 ...
         -0.004364327211358 -0.010013251577232 -0.021344241245400 -0.042182820580118 -0.077080889653194...
         -0.129740392318591 -0.200064921294891 -0.280328573340852 -0.352139052257134 -0.386867664739069 ...
         -0.351974030208595 -0.223363323458050 0 0.286427448595213 0.574058766243311 ...
         0.788100265785590 0.867325070584078 0.788100265785590 0.574058766243311 0.286427448595213 0 ...
         -0.223363323458050 -0.351974030208595 -0.386867664739069 -0.352139052257134...
         -0.280328573340852 -0.200064921294891 -0.129740392318591 -0.077080889653194 -0.042182820580118 ...
         -0.021344241245400 -0.010013251577232 -0.004364327211358 -0.001770119249103 -6.689305101192819e-04...
         -2.357742589814283e-04 -7.757327341237223e-05];

    b1 = resample(b1,fs,250);
    bpfecg = filtfilt(b1,1,ecg)';
    
    if (sum(abs(ecg-median(ecg))>MIN_AMP)/NB_SAMP)>0.05
        % if 20% of the samples have an absolute amplitude which is higher
        % than MIN_AMP then we are good to go.
        
        % == P&T operations
        dffecg = diff(bpfecg');  % (4) differentiate (one datum shorter)
        sqrecg = dffecg.*dffecg; % (5) square ecg
        intecg = filter(ones(1,INT_NB_COEFF),1,sqrecg); % (6) integrate
        mdfint = medfilt1(intecg,MED_SMOOTH_NB_COEFF);  % (7) smooth
        delay  = ceil(INT_NB_COEFF/2); 
        mdfint = circshift(mdfint,-delay); % remove filter delay for scanning back through ECG

        % look for some measure of signal quality with signal fid_vec? (FIXME)
        if isempty(fid_vec); mdfintFidel = mdfint; else mdfintFidel(fid_vec>2) = 0; end;

        % == P&T threshold
        if NB_SAMP/fs>90; xs=sort(mdfintFidel(fs:fs*90)); else xs = sort(mdfintFidel(fs:end)); end;

        if isempty(MAX_FORCE)
           if NB_SAMP/fs>10
                ind_xs = ceil(98/100*length(xs)); 
                en_thres = xs(ind_xs); % if more than ten seconds of ecg then 98% CI
            else
                ind_xs = ceil(99/100*length(xs)); 
                en_thres = xs(ind_xs); % else 99% CI  
            end 
        else
           en_thres = MAX_FORCE;
        end

        % build an array of segments to look into
        poss_reg = mdfint>(THRES*en_thres); 

        % in case empty because force threshold and crap in the signal
        if isempty(poss_reg); poss_reg(10) = 1; end;

        % == P&T QRS detection & search back
        if SEARCH_BACK
            indAboveThreshold = find(poss_reg); % ind of samples above threshold
            RRv = diff(tm(indAboveThreshold));  % compute RRv
            medRRv = median(RRv(RRv>0.01));
            indMissedBeat = find(RRv>1.5*medRRv); % missed a peak?
            % find interval onto which a beat might have been missed
            indStart = indAboveThreshold(indMissedBeat);
            indEnd = indAboveThreshold(indMissedBeat+1);

            for i=1:length(indStart)
                % look for a peak on this interval by lowering the energy threshold
                poss_reg(indStart(i):indEnd(i)) = mdfint(indStart(i):indEnd(i))>(0.5*THRES*en_thres);
            end
        end

        % find indices into boudaries of each segment
        left  = find(diff([0 poss_reg'])==1);  % remember to zero pad at start
        right = find(diff([poss_reg' 0])==-1); % remember to zero pad at end

        % looking for max/min?
        if SIGN_FORCE
            sign = SIGN_FORCE;
        else
            nb_s = length(left<30*fs);
            loc  = zeros(1,nb_s);
            for j=1:nb_s
                [~,loc(j)] = max(abs(bpfecg(left(j):right(j))));
                loc(j) = loc(j)-1+left(j);
            end
            sign = mean(ecg(loc));  % FIXME: change to median?  
        end

        % loop through all possibilities 
        compt=1;
        NB_PEAKS = length(left);
        maxval = zeros(1,NB_PEAKS);
        maxloc = zeros(1,NB_PEAKS);
        for i=1:NB_PEAKS
            if sign>0
                % if sign is positive then look for positive peaks
                [maxval(compt) maxloc(compt)] = max(ecg(left(i):right(i)));
            else
                % if sign is negative then look for negative peaks
                [maxval(compt) maxloc(compt)] = min(ecg(left(i):right(i)));
            end
            maxloc(compt) = maxloc(compt)-1+left(i); % add offset of present location

            % refractory period - has proved to improve results
            if compt>1
                if maxloc(compt)-maxloc(compt-1)<fs*REF_PERIOD && abs(maxval(compt))<abs(maxval(compt-1))
                    maxloc(compt)=[]; maxval(compt)=[];
                elseif maxloc(compt)-maxloc(compt-1)<fs*REF_PERIOD && abs(maxval(compt))>=abs(maxval(compt-1))
                    maxloc(compt-1)=[]; maxval(compt-1)=[];
                else
                    compt=compt+1;
                end
            else
                % if first peak then increment
                compt=compt+1;
            end
        end

        qrs_pos = maxloc; % datapoints QRS positions 
        R_t = tm(maxloc); % timestamps QRS positions
        R_amp = maxval; % amplitude at QRS positions
        hrv = 60./diff(R_t); % heart rate
    else
        % this is a flat line
        qrs_pos = [];
        R_t = [];
        R_amp = [];
        hrv = [];
        sign = [];
        en_thres = [];
    end
catch ME
    rethrow(ME);
    for enb=1:length(ME.stack); disp(ME.stack(enb)); end;
    qrs_pos = [1 10 20]; sign = 1; en_thres = 0.5; 
end

% == plots
if debug
    figure;
    FONTSIZE = 20;
    ax(1) = subplot(4,1,1); plot(tm,ecg); hold on;plot(tm,bpfecg,'r')
        title('raw ECG (blue) and zero-pahse FIR filtered ECG (red)'); ylabel('ECG');
        xlim([0 tm(end)]);  hold off;
    ax(2) = subplot(4,1,2); plot(tm(1:length(mdfint)),mdfint);hold on;
        plot(tm,max(mdfint)*bpfecg/(2*max(bpfecg)),'r',tm(left),mdfint(left),'og',tm(right),mdfint(right),'om'); 
        title('Integrated ecg with scan boundaries over scaled ECG');
        ylabel('Int ECG'); xlim([0 tm(end)]); hold off;
    ax(3) = subplot(4,1,3); plot(tm,bpfecg,'r');hold on;
        plot(R_t,R_amp,'+k');
        title('ECG with R-peaks (black) and S-points (green) over ECG')
        ylabel('ECG+R+S'); xlim([0 tm(end)]); hold off;
    ax(4) = subplot(4,1,4); plot(R_t(1:length(hrv)),hrv,'r+')
        hold on, title('HR')
        ylabel('RR (s)'); xlim([0 tm(end)]);
    
    %linkaxes(ax,'x');
    set(gca,'FontSize',FONTSIZE);
    allAxesInFigure = findall(gcf,'type','axes');
    set(allAxesInFigure,'fontSize',FONTSIZE);
end


% NOTES
%   Finding the P&T energy threshold: in order to avoid crash due to local 
%   huge bumps, threshold is choosen at 98-99% of amplitude distribution. 
%   first sec removed for choosing the thres because of filter init lag.
%   
%   Search back: look for missed peaks by lowering the threshold in area where the 
%   RR interval variability (RRv) is higher than 1.5*medianRRv
% 
%   Sign of the QRS (signForce): look for the mean sign of the R-peak over the
%   first 30sec when looking for max of abs value. Then look for the
%   R-peaks over the whole record that have this given sign. This allows to
%   not alternate between positive and negative detections which might
%   happen in some occasion depending on the ECG morphology. It is also
%   better than forcing to look for a max or min systematically.











 end

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

function en=Teager(sig)


Lsig  = length(sig);

en = sig(2:Lsig-1).^2 - sig(1:Lsig-2).*sig(3:Lsig);

end

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

function [e,se,A,B]=sampen(y,M,r,sflag,cflag,vflag)
%function e=sampen(y,M,r);
%
%Input Parameters
%
%y input signal vector
%M maximum template length (default M=5)
%r matching threshold (default r=.2)
%
%Output Parameters
%
%e sample entropy estimates for m=0,1,...,M-1
%
%Full usage:
%
%[e,se,A,B]=sampen(y,m,r,sflag,cflag,vflag)
%
%Input Parameters
%
%sflag flag to standardize signal(default yes/sflag=1)
%cflag flag to use fast C code (default yes/cflag=1)
%vflag flag to calculate standard errors (default no/vflag=0)
%
%Output Parameters
%
%se standard error estimates for m=0,1,...,M-1
%A number of matches for m=1,...,M
%B number of matches for m=0,...,M-1
% (excluding last point in Matlab version)
if ~exist('M')|isempty(M),M=5;end
if ~exist('r')|isempty(r),r=.2;end
if ~exist('sflag')|isempty(sflag),sflag=1;end
if ~exist('cflag')|isempty(cflag),cflag=1;end
if ~exist('vflag')|isempty(cflag),vflag=0;end
y=y(:);
n=length(y);
if sflag>0
y=y-mean(y);
s=sqrt(mean(y.^2));
y=y/s;
end
if nargout>1
if vflag>0
se=sampense(y,M,r);
else
se=[];
end
end
if cflag>0
[match,R]=cmatches(y,n,r);
match=double(match);
else
[e,A,B]=sampenc(y,M,r);
return
end
k=length(match);
if k<M
match((k+1):M)=0;
end
N=n*(n-1)/2;
A=match(1:M);
B=[N;A(1:(M-1))];
N=n*(n-1)/2;
p=A./B;
e=-log(p);
end

% rePE.m - algorithm for the fast calculation of a robust empirical
% permutaion entropy in maximally overlapping sliding windows

% INPUT (x - the considered time series, Tau - a delay, 
% d - an order of the ordinal patterns, WS - size of a sliding window,
% thr1 and thr2 - the lower and upper thresholds)
% OUTPUT [ re_PE - the values of robust empirical permutation entropy,
% MD - the values of MD]

function [MD, re_PE] = rePE(x, Tau, d, WS, thr1, thr2)
load(['table' num2str(d) '.mat']); % the precomputed table
pTbl = eval(['table' num2str(d)]);    
Length = numel(x);                 % length of the time series
d1 = d+1; 
dTau = d*Tau; 
nPat = factorial(d1);              % amount of ordinal patterns of order d               
opd = zeros(1, nPat);              % distribution of ordinal patterns
op = zeros(1, d);                  % ordinal pattern $(i_1,i_2,...,i_d)$
ancNum = nPat./factorial(2:d1);    % ancillary numbers       
MDthr = (d+1)*d/8;
prevOP = zeros(1, Tau);            % previous ordinal patterns for $1:\tau$
opW = zeros(1, WS);                % ordinal patterns in the window
re_PE = zeros(1, Length- WS-dTau);

MD = zeros(1, Length);
for iTau = 1:Tau
    MDar1 = zeros(1, d);
    MDar2 = zeros(1, d);
    for i  = 1:d
        MDar1(i)=sum(abs(x(iTau+(i-1)*Tau)-x(iTau+i*Tau:Tau:iTau+dTau))<thr1);
        MDar2(i)=sum(abs(x(iTau+(i-1)*Tau)-x(iTau+i*Tau:Tau:iTau+dTau))>thr2);
    end
    MD(iTau) = sum(MDar1)+sum(MDar2);
    MDar1(1:d-1) = MDar1(2:d);
    MDar2(1:d-1) = MDar2(2:d);
    MDar1(d) = 0;
    MDar2(d) = 0;
    for i = iTau+Tau:Tau:Length-dTau-Tau
        for j =0:d-1
            MDar1(j+1) = MDar1(j+1) + (abs( x(i+j*Tau)-x(i+dTau) ) < thr1);
            MDar2(j+1) = MDar2(j+1) + (abs( x(i+j*Tau)-x(i+dTau) ) > thr2);
        end
        MD(i) = sum(MDar1)+sum(MDar2);
        MDar1(1:d-1) = MDar1(2:d);
        MDar1(d) = 0;
        MDar2(1:d-1) = MDar2(2:d);
        MDar2(d) = 0;
    end
end

for iTau = 1:Tau                     % the first sliding window
  cnt = iTau; 
  op(1) = (x(dTau+iTau-Tau) >= x(dTau+iTau));
  for j = 2:d
    op(j) = sum(x((d-j)*Tau+iTau) >= x((d1-j)*Tau+iTau:Tau:dTau+iTau));
  end        
  opW(cnt) = sum(op.*ancNum);        % the first ordinal pattern
  OPnumber = opW(cnt);              
  if(MD(cnt)<MDthr)
        opd(OPnumber+1) = opd(OPnumber+1)+1;  
  end
  for j = dTau+Tau+iTau:Tau:WS+dTau  % loop for the first window
    cnt = cnt+Tau;                               
    posL = 1;                        % the position $l$ of the next point
    for i = j-dTau:Tau:j-Tau
        if(x(i) >= x(j)) 
          posL = posL+1; 
        end
    end  
    opW(cnt) = pTbl(opW(cnt-Tau)*d1+posL);
    OPnumber = opW(cnt);
    if(MD(cnt)<MDthr)
        opd(OPnumber+1) = opd(OPnumber+1)+1;            
    end
  end 
  prevOP(iTau) = opW(cnt);
end    
ordDistNorm = opd/sum(opd);
re_PE(WS+Tau*d) = -nansum(ordDistNorm(1:nPat).*log(ordDistNorm(1:nPat)))/d;       

iTau = mod(WS, Tau)+1;          % current shift $1:\tau$
iPat = 1;                       % position of the current pattern in the window
for t = WS+Tau*d+1:Length       % loop over all points
  posL = 1;                     % the position $l$ of the next point
  for j = t-dTau:Tau:t-Tau
    if(x(j) >= x(t)) 
        posL = posL+1; 
    end
  end                          
  nNew = pTbl(prevOP(iTau)*d1+posL); % "incoming" ordinal pattern         
  nOut = opW(iPat);                  % "outcoming" ordinal pattern 
  prevOP(iTau) = nNew;
  opW(iPat) = nNew; 
  nNew = nNew+1;
  nOut = nOut+1;       
  % update the distribution of ordinal patterns 
  if (MD(t-dTau)<MDthr)   
     opd(nNew) = opd(nNew)+1; % "incoming" ordinal pattern
  end
  if (MD(t-WS-dTau)<MDthr)
     opd(nOut) = opd(nOut)-1; % "outcoming" ordinal pattern
  end
ordDistNorm = opd/sum(opd);
re_PE(t) = -nansum(ordDistNorm(1:nPat).*log(ordDistNorm(1:nPat)))/d;     
  
  iTau = iTau+1; 
  iPat = iPat+1;
  if(iTau > Tau) iTau = 1; end
  if(iPat > WS) iPat = 1; end
end 
re_PE = re_PE(WS+Tau*d:end);
end

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