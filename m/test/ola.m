close all
clear all
h = [1,7,8];
N = length(h);

x = [5,8,9,6,3,4,8,2,1,7,5,6];
M = 4;

L = M+N-1;


x1 = x(1:4);
x2 = x(5:8);
x3 = x(9:12);

H = fft(h,L);

Y1 = fft(x1,L).*H;
Y2 = fft(x2,L).*H;
Y3 = fft(x3,L).*H;

y1 = ifft(Y1,L);
y2 = ifft(Y2,L);
y3 = ifft(Y3,L);

y_conv = conv(h,x)
y_ola = [y1(1:4),y1(5:6)+...
                y2(1:2),y2(3:4),y2(5:6)+...
                                 y3(1:2),y3(3:6)]
stem(y_conv);
hold on;
stem(y_ola);