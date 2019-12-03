close all;

total_input_f=aec('rec_mic_0_0.pcm', 'rec_mic_1_0.pcm', 'rec_spk_l_0.pcm');
tmp1 = total_input_f(:,1:129);
tmp2 = total_input_f(:,130:258);
[out1, out2] = compose(tmp1, tmp2);
x = size(out1, 1);
y = size(out1, 2);
out1 = reshape(out1, [1, x*y]);
out2 = reshape(out2, [1, x*y]);

file_id=fopen('out_aec_0.pcm','wb');
fwrite(file_id, out1,'int16');
fclose(file_id);

file_id=fopen('out_aec_1.pcm','wb');
fwrite(file_id, out2,'int16');
fclose(file_id);
