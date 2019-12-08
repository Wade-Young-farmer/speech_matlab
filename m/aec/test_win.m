clc;
close all;
clear all;
% 这不是sin函数求幅频响应的方式
t=-5:0.01:4.99;
x=sin(2*pi*t);
plot(t, x);
X=fft(x);
Y=fft(x(1:50));
Z=fft(x(1:100));
figure;
plot(abs(X));
hold on;
plot(abs(Y), '*');
hold on;
plot(abs(Z));

% 这不是窗函数的幅频响应求的方式
x1=zeros(1, 100);
x1(40:1:59)=1;
y1=fft(x1);
figure;
plot(abs(y1));

%{
--------------------------------------------------------------------------- 
File:Matlab的窗函数,矩形窗                               
功能：降低旁瓣水平
参数： 
--------------------------------------------------------------------------- 
%}
%N =51 
%==========================================================================
%求矩形窗的频率响应图  
%==========================================================================
W = linspace(-pi,pi,4096);
wn0 = rectwin(20);   %矩形窗函数 
%20*log10(abs(WN))  
[h1,w0] = freqz(wn0,1,W); 
%subplotfigure(5,1,1);  
figure;
plot(w0/pi,20*log10(abs(h1/max(h1))));  
axis([-1 1 -100 0]); 
xlabel('归一化频率 /\pi');  
ylabel('20log_{10}|W(e^{j\omega})| /dB');  
title('矩形窗的傅里叶变换'); 
set(gca,'YTick',[-100 -80 -60 -40 -20 0])  
set(gca,'XTick',[-1 :0.2: 1])   
%set(gca,'XAxisLocation','top');%设置X轴在上方  
set(gca,'XAxisLocation','bottom');%设置X轴在下方  
set(gca,'YAxisLocation','left'); %设置Y轴在左方  
text(1,-124,'\pi');%gtext('\pi');
figure;
plot(w0/pi, phase(h1));
