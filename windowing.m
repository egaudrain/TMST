function windows = windowing(signal, sample_rate, width, shift, gwin)
% windows = WINDOWING(signal, sample_rate, width, overlap, gwin) applies a
% windowing (with parameters width and shift) to an input signal with
% sample_rate. Shift is the shift between successive windows, in sec. If
% gwin is 1, the function uses a Gaussian window.
%
% Leo Varnet - 08/2023

signal = signal(:)';
width_in_sample = floor(width*sample_rate);
shift_in_sample = floor(shift*sample_rate);

%windowing function
if gwin
    G = gausswin(width_in_sample);
    norm = 1;%/(((width_in_sample-1)/(2*2.5))*sqrt(2*pi));
else
    G = 1;%ones(1,width_in_sample);
    norm = 1;
end

i=1;
windows = [];
while i+width_in_sample<length(signal)
    windows = [windows, (signal(i:i+width_in_sample-1).*G')'/norm];
    i=i+shift_in_sample;
end
end

