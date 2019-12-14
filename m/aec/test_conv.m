close all;
clear all;
clc;
t_size=128;
t=0:1/t_size:1-1/t_size;
% h=sin(pi*t);
% h=ones(size(t));
h=hamming(t_size);
delta=1;
unit=[1, 1];
h1 = conv(h,delta);
h2 = conv(h,unit);

delta_w = fft(delta);
h_w =fft(h);
H_W = h_w .* delta_w;

W = 0:2*pi/t_size:2*pi*(t_size-1)/t_size; 
H_W1 = freqz(h, 1, W);

W2=linspace(0, 2*pi, 2* t_size + 1);
H_W2 = freqz(h, 1, W2);
% 证明了窗函数的fft 等价于
% 直接求幅频响应，只是画幅频响应的时候，w的取值会比较精细，如果直接按fft的点数来取值的话，两者的曲线就完全吻合
figure;
plot(W/(2*pi), 20*log10(abs(H_W)/max(abs(H_W))), 'b');
hold on;
plot(W/(2*pi), 20*log10(abs(H_W1)/max(abs(H_W1))), 'g');
hold on;
plot(W2/(2*pi), 20*log10(abs(H_W2)/max(abs(H_W2))), 'r');

%%%%%%%%%%%%%%%%%
t_size=128;
t=0:1/t_size:2-1/t_size;
% h=sin(pi*t);
% h=ones(size(t));
h=hamming(t_size * 2);
x = 2* sin(2*pi*t) + sin(3*pi*t +0.5);
x_w = fft(x);
x_h = x .* h.';
X_H_W = fft(x_h);

% 证明了， 此刻t的长度有一个周期，不需要加窗，加窗反而会导致频谱泄露
figure;
plot(abs(x_w), '*');
hold on;
stem(abs(X_H_W));

%%%%%%%%%%%%%%%%%
t_size=128;
t=0:1/t_size:3-1/t_size;
% h=sin(pi*t);
% h=ones(size(t));
h=hamming(t_size * 3);
x = 2* sin(2*pi*t) + sin(3*pi*t +0.5);
x_w = fft(x);
x_h = x .* h.';
X_H_W = fft(x_h);

% 证明了，此刻t的长度有一个半周期，会导致频谱泄露，加窗可以缓解频谱泄露
figure;
plot(abs(x_w), '*');
hold on;
stem(abs(X_H_W));

X_H_W1 = 1 / ( 2* pi) * conv(x_w, h_w);
% 证明了信号和窗函数分别做完fft后，再做频域卷积，再除以2pi, 不等价于
% 加窗后再做fft, 首先卷积出来的频谱长度就和目标fft的长度不相等， 如果fft做完zero padding之后可能会相等？
figure;
plot(abs(X_H_W1));
hold on;
plot(abs(X_H_W), 'g');

