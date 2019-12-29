close all
clear all
h = [1,7,8];
N = length(h);

x0 = [5,8,9,6,3,4,8,2,1,7,5,6];
M = 4;

x = [zeros(1,N-1),x0];

L = M+N-1;


x1 = x(1:6);
x2 = x(5:10);
x3 = x(9:14);
x4 = [x(13:14),zeros(1,4)];

H = fft(h,L);

Y1 = fft(x1,L).*H;
Y2 = fft(x2,L).*H;
Y3 = fft(x3,L).*H;
Y4 = fft(x4,L).*H;


y1 = ifft(Y1,L);
y2 = ifft(Y2,L);
y3 = ifft(Y3,L);
y4 = ifft(Y4,L);


y_conv = conv(h,x0)
y_ols = [y1(3:6),y2(3:6),y3(3:6),y4(3:6)]
stem(y_conv);
hold on;
stem(y_ols);


