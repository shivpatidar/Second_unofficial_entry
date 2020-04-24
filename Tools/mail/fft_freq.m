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
