clear all;
clc;
close all;

run para
f=fopen(fle,'r');
x=fread(f,inf,'int16');
fclose(f);
figure;
plot(x);
x=hpf(x,fs);
% sound(x, fs)
out_put=subband_analyze(x,fs,real_win);
figure;
plot(abs(out_put(:,3)));
test=out_put(:,3).';
