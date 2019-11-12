function my_agc(filename,vad_flag,outfile)
Srate=16000;
file_id=fopen(filename, 'r');
x=fread(file_id, inf, 'int16');
fclose(file_id);

enframed_x=enframe(x,128);
gain_buffer=zeros(1,10);
gain_curr = 1.0;
gain_pre = 1.0;
gain_last = 1.0;
output_buf=zeros(1,64);
for i=1:size(enframed_x,1)
    input = enframed_x(i,:);
    peak = max(input);
    energy = sqrt(mean(input .* input));
    gain_1 = 16218.83752939778 /(energy + 1);
    gain_2 = 26000/(peak+1);
    gain_3 = 24000/(peak+1); % 唤醒路此处为26000
    gain = min([gain_1, gain_3, gain_2/1.5]);
    gain_clip = gain_2;
    narrow_factor = 1.0;
    gain_buffer(mod(i-1, 10)+1)=gain*narrow_factor;
    %% audio buffer smooth
    tmp=input(65:128);
    input(65:128)=input(1:64);
    input(1:64)=output_buf;
    output_buf=tmp;
    if i > 10 && vad_flag == 1
        gain_tmp=sort(gain_buffer);
        gain_ = gain_tmp(2); %找一个次小值
        % gain_ smooth 过程, gain_, gain_last, gain_final的值应该经过平滑
        gain_final = gain_;
        
        gain_curr = min(gain_final, gain_clip);
        gain_pre = min(gain_pre, gain_clip);
        for j = 1:128
            input(j) = input(j) * (gain_pre + (gain_curr-gain_pre)*(j-1)/128);
        end
        gain_pre = gain_pre + 127 /128 * (gain_curr-gain_pre);
        gain_last = gain_pre;
    else
        gain_curr = min(gain_last, gain_clip);
        gain_pre = min(gain_pre, gain_clip);
        for j = 1:128
            input(j) = input(j) * (gain_pre + (gain_curr-gain_pre)*(j-1)/128);
        end
        gain_pre = gain_curr;
    end
    
    xfinal((i-1)*128+1:i*128)=input;
end

% xfinal=x;
file_id=fopen(outfile,'wb');
fwrite(file_id,xfinal,'int16');
fclose(file_id);
end