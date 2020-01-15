clc
clear all

N=128;
% The parameters of the first filter
%sos1 = [1 -2 1 1 -1.95998 0.9615];
%gain1 = 0.0013094*4096;

sos1 = [1 -2 1 1 -1.95998 0.9615];
gain1 = 1.2512;

% The parameters of the second filter
sos2 = [2 -2 0 1 -0.96148 0];
gain2 = 367.1438/4096;

sos_final = [sos1;sos2];
gain_final = gain1 * gain2;
[b, a] = sos2tf(sos_final, gain_final/8);
equal_trans = [b,a];

[H1,w] = freqz(sos1*gain1,N);
figure(1);
subplot(2,1,1);
plot(w/pi, 20*log10(abs(H1)));grid;
%plot(w/pi,abs(H1));grid;
title('Frequency response of the first filter');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Magnitude (dB)');

subplot(2,1,2);
plot(w/pi, angle(H1));grid;
title('Frequency response of the first filter');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Phase');


[H2,w] = freqz(sos2*gain2,N);
figure(2);
subplot(2,1,1);
plot(w/pi,20*log10(abs(H2)));grid;
%plot(w/pi,abs(H2));grid;
title('Frequency response of the second filter');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Magnitude (dB)');

subplot(2,1,2);
plot(w/pi,angle(H2));grid;
title('Frequency response of the second filter');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Phase');

[H3,w] = freqz(equal_trans,N);
figure(3);
subplot(2,1,1);
plot(w/pi,20*log10(abs(H3)));grid;
%plot(w/pi,abs(H3));grid;
title('Frequency response of the final filter');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Magnitude (dB)');
subplot(2,1,2);
plot(w/pi,angle(H3));grid;
title('Frequency response of the final filter');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Phase');

