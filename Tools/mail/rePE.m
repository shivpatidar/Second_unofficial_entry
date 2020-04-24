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