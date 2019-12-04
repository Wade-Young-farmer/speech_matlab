function [x_enframe, x2_enframe, r_enframe, x_f, x2_f, r_f]=hpf(mic_name, mic_name_2, ref_name, filter_coeff)
file_id=fopen(mic_name, 'r');
x=fread(file_id, inf, 'int16');
fclose(file_id);

file_id=fopen(mic_name_2, 'r');
x2=fread(file_id, inf, 'int16');
fclose(file_id);

file_id=fopen(ref_name, 'r');
r=fread(file_id, inf, 'int16');
fclose(file_id);

x_enframe = enframe(x,128);
x2_enframe = enframe(x2,128);
r_enframe = enframe(r,128);

frame_size=size(x_enframe, 1);
x_f=zeros(frame_size, 129);
x2_f=zeros(frame_size, 129);
r_f=zeros(frame_size, 129);

b=[1, -2, 1];
a=[1, -1.95998, 0.9615];
d=[2, -2];
c=[1, -0.96148];
g = [0.0013094, 367.1433];
b0 = g(1) * g(2) * conv(b, d);
a0 = conv(a, c);

buf_1 = zeros(1, 768);
buf_2 = zeros(1, 768);
buf_r = zeros(1, 768);

for i=1:frame_size
    input = x_enframe(i,:);
    input = filter(b0, a0, input);
    x_enframe(i,:)=input;
    
    input2=x2_enframe(i,:);
    input2=filter(b0, a0, input2);
    x2_enframe(i,:)=input2;
    
    ref=r_enframe(i,:);
    ref=filter(b0, a0, ref);
    r_enframe(i,:)=ref;
    
    buf_1(641:768)=input;
    temp = filter_coeff .* buf_1;
    temp2 = zeros(1, 256);
    for j = 1:256
        temp2(j) = temp(j) + temp(j + 256) + temp(j + 512);
    end
    buf_1(1:640)=buf_1(129:768);
    input_f=fft(temp2);
    input_f=input_f(1:129);
    x_f(i,:)=input_f;
    
    buf_2(641:768)=input2;
    temp = filter_coeff .* buf_2;
    temp2 = zeros(1, 256);
    for j = 1:256
        temp2(j) = temp(j) + temp(j + 256) + temp(j + 512);
    end
    buf_2(1:640)=buf_2(129:768);
    input_f=fft(temp2);
    input_f=input_f(1:129);
    x2_f(i,:)=input_f;
    
    buf_r(641:768)=ref;
    temp = filter_coeff .* buf_r;
    temp2 = zeros(1, 256);
    for j = 1:256
        temp2(j) = temp(j) + temp(j + 256) + temp(j + 512);
    end
    buf_r(1:640)=buf_r(129:768);
    input_f=fft(temp2);
    input_f=input_f(1:129);
    r_f(i,:)=input_f;
end

end