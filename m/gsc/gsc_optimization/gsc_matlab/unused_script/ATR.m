
function H=ATR(v1,v2,v3,v4,v5,v6)   %输入为6路语音对应频点数
frm_l=1024;          %帧长
frm_lap=512;        %帧移
use_fnum=20;
voice1=enframe(v1,hanning(frm_l,'periodic'),frm_lap);         %voice1~6分别存6路语音，为矩阵，每一行对应一帧语音
voice2=enframe(v2,hanning(frm_l,'periodic'),frm_lap);  
voice3=enframe(v3,hanning(frm_l,'periodic'),frm_lap);  
voice4=enframe(v4,hanning(frm_l,'periodic'),frm_lap);  
voice5=enframe(v5,hanning(frm_l,'periodic'),frm_lap);  
voice6=enframe(v6,hanning(frm_l,'periodic'),frm_lap); 

fai1=zeros(6,frm_l/2+1);      %变量初始化
fai2=zeros(1,frm_l/2+1); 
fai3=zeros(6,frm_l/2+1); 
fai4=zeros(1,frm_l/2+1); 

for i=1:use_fnum   %控制帧循环
    
        fv1=fft(voice1(i,:));        %进行傅里叶变换，补零到1024点
        fv2=fft(voice2(i,:));
        fv3=fft(voice3(i,:));
        fv4=fft(voice4(i,:));
        fv5=fft(voice5(i,:));
        fv6=fft(voice6(i,:));
        fv=[fv1;fv2;fv3;fv4;fv5;fv6];
    
    fai2=fai2+fv(1,1:frm_l/2+1).*conj(fv(1,1:frm_l/2+1));
    fai4=fai4+fv(1,1:frm_l/2+1).*conj(fv(1,1:frm_l/2+1)).*fv(1,1:frm_l/2+1).*conj(fv(1,1:frm_l/2+1));
    
    for m=1:6
        fai1(m,:)=fai1(m,:)+fv(1,1:frm_l/2+1).*conj(fv(1,1:frm_l/2+1)).*fv(m,1:frm_l/2+1).*conj(fv(1,1:frm_l/2+1));
        fai3(m,:)=fai3(m,:)+fv(m,1:frm_l/2+1).*conj(fv(1,1:frm_l/2+1));
    end
end
fai1=fai1/use_fnum;
fai2=fai2/use_fnum;
fai3=fai3/use_fnum;
fai4=fai4/use_fnum;
fai5=fai2.*conj(fai2);
for m=1:6
    H(m,:)=(fai1(m,:)-fai2.*fai3(m,:))./(fai4-fai5);
end
end