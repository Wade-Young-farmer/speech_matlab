function output=subband_analyze(signal,fs,win)
addpath('../basic_tbx');
time=0.008;
frame_len=time*fs;
inc=frame_len;
fake_win=ones(1, frame_len*6);
assert(length(win)==(frame_len*6));
x_temp=enframe(signal,fake_win,inc).';
x_c=zeros(frame_len*6,5);
for i=1:5
    x_c_temp=[zeros(frame_len*(6-i),1);x_temp(1:frame_len*i,1)]; 
    x_c_temp=x_c_temp .* (win.');
    x_c(:,i)=x_c_temp;
end

x=enframe(signal,win,inc).';
x=[x_c,x];
frame_num=size(x, 2);
xx=zeros(frame_len*2,frame_num);
for i=1:frame_num
    %Ê±ÓòµÄµþ¼Ó
    real_x=x(1:frame_len*2,i)+x(frame_len*2+1:frame_len*4,i)+x(frame_len*4+1:frame_len*6,i);
    xx(:,i)=real_x;
end
y=fft(xx);
output=y(1:frame_len,:);




