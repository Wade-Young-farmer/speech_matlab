function [x_out, x2_out]=compose(input1, input2, filter_coeff)
frame_size=size(input1, 1);
x_out=zeros(frame_size, 128);
x2_out=zeros(frame_size, 128);

buf_1=zeros(1, 768);
buf_2=zeros(1, 768);

for i=1:frame_size
    f1=input1(i,:);
    temp=zeros(1,256);
    temp(1:129)=f1;
    temp(130:256)=conj(f1(128:-1:2));
    t1=ifft(temp) * 256;
    temp2 = [t1, t1, t1];
    buf_1 = buf_1 + filter_coeff .* temp2;
    output1 = buf_1(1:128)*128;
    x_out(i, :)=output1;
    buf_1(1:640)=buf_1(129:768);
    buf_1(641:768)=0;
    
    f2=input2(i,:);
    temp=zeros(1,256);
    temp(1:129)=f2;
    temp(130:256)=conj(f2(128:-1:2));
    t2=ifft(temp) * 256;
    temp2 = [t2, t2, t2];
    buf_2 = buf_2 + filter_coeff .* temp2;
    output2 = buf_2(1:128)*128;
    x2_out(i, :)=output2;
    buf_2(1:640)=buf_2(129:768);
    buf_2(641:768)=0;    
end
end