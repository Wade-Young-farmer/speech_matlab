clear all;
clc;
close all;

run para
fle='data/rec_spk_l.pcm';
f=fopen(fle,'r');
x=fread(f,inf,'int16');
fclose(f);
%figure;
%plot(x);
fle='data/rec_mic_0.pcm';
f=fopen(fle, 'r');
x1=fread(f,inf,'int16');
fclose(f);
fle='data/rec_mic_1.pcm';
f=fopen(fle, 'r');
x2=fread(f,inf,'int16');
fclose(f);

x=hpf(x,fs);
% sound(x, fs)
xr=subband_analyze(x,fs,real_win);

x1=hpf(x1,fs);
x1=subband_analyze(x1,fs,real_win);

x2=hpf(x2,fs);
x2=subband_analyze(x2,fs,real_win);

figure;
plot(abs(x1(:,1)));
test=x1(:,1).';
