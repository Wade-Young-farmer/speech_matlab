close all;
clc;
x1 = [1,7,8,9,5,4,6,3,2];
x2 = [5,8,9,6,3,4,8,2,1,7,5,6,7];
y = conv(x1,x2);
size(y, 2)
x3 = [x1, zeros(1, size(y,2)-size(x1, 2))];
x4 = [x2, zeros(1, size(y,2)-size(x2, 2))];
y1 = cconv(x3, x4, size(y, 2));

stem(y)
hold on;
stem(y1, 'g');
% L = 15;
% hold on;
% stem(ifft(fft(x1,L).*fft(x2,L)))