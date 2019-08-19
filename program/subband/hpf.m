%% 每帧都做高通滤波， iir滤波器，零点和极点配置，类似于图9.3.1直接I型实现
function output=hpf(signal,fs)
addpath('../basic_tbx');
b=[1.0, -2.0, 1.0];
a=[1, -1.95998,  0.9615];
d=[2, -2];
c=[1.0, -0.96148];
g=[0.0013094, 367.1433];

time=0.008;
frame_len=time*fs;
inc=frame_len;
fake_win=ones(1,frame_len);

x_temp=enframe(signal,fake_win,inc).';
frame_num=size(x_temp,2);
x_temp=filter(b,a,x_temp);
x_temp=x_temp * g(1);
x_temp=filter(d,c,x_temp);
x_temp=x_temp*g(2);
output=reshape(x_temp,[1, frame_num*frame_len]);

