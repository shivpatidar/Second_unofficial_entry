function [filt_signal1] = BP_filter_ECG(ecg,fs)

ecg=ecg;
fs=fs;

d = designfilt('bandpassiir','FilterOrder',6, ...
    'HalfPowerFrequency1',1,'HalfPowerFrequency2',35, ...
    'SampleRate',fs);
            

% Notch filter
%wo = 60/(fs/2);  bw = wo/2;
%[b,a] = iirnotch(wo,bw);
%y = filter(b,a,ecgn);
                                % plot(ecgn)
                                % hold on
                                % plot(y)

                                % plot(abs(fft(x1)))
                                % hold on
                                % plot(abs(fft(y)))
                                % 
                                %fvtool(b,a);

%% Filtering 

    %yy_1 = filter(b,a,xx_1);
    %yy_2=filtfilt(d,yy_1);
    
    filt_signal1=filtfilt(d,ecg);
   
end

