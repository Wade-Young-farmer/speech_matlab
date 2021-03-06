%
% pr6_2_3 
clear all; clc; close all;

filedir=[];                             % 指定文件路径
filename='bluesky1.wav';                % 指定文件名
fle=[filedir filename]                  % 构成路径和文件名的字符串
[xx,fs]=wavread(fle);                   % 读入数据文件
x=xx/max(abs(xx));                      % 幅度归一化
N=length(xx);                           % 取信号长度
time=(0:N-1)/fs;                        % 计算时间刻度

wlen=200; inc=80;                       % 设置帧长和帧移
IS=0.25; overlap=wlen-inc;              % 设置前导无话段长度
NIS=fix((IS*fs-wlen)/inc +1);           % 计算前导无话段帧数
y=enframe(x,wlen,inc)';                 % 分帧
etemp=sum(y.^2);                        % 求取短时平均能量
etemp=etemp/max(etemp);                 % 能量幅值归一化
fn=size(y,2);                           % 帧数
T1=0.002;                               % 设置阈值
T2=0.01;
frameTime=frame2time(fn, wlen, inc, fs);% 计算各帧对应的时间
[voiceseg,vsl,SF,NF]=vad_param1D(etemp,T1,T2);% 用一个参数端点检测
% 作图
subplot 211; plot(time,x,'k'); hold on
title('纯语音男声“蓝天，白云，碧绿的大海”波形');
ylabel('幅值'); axis([0 max(time) -1 1]); 
for k=1 : vsl
    nx1=voiceseg(k).begin; nx2=voiceseg(k).end;
    fprintf('%4d   %4d   %4d\n',k,nx1,nx2);
    line([frameTime(nx1) frameTime(nx1)],[-1 1],'color','k','LineStyle','-');
    line([frameTime(nx2) frameTime(nx2)],[-1 1],'color','k','LineStyle','--');
end
subplot 212; plot(frameTime,etemp,'k');
title('语音短时能量图');
ylabel('幅值'); axis([0 max(time) 0 1]);
xlabel('时间/s');
line([0 max(time)],[T1 T1],'color','k','LineStyle','-');
line([0 max(time)],[T2 T2],'color','k','LineStyle','--');
