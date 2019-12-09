close all;
clear all;
token = [4.082482904638630e-001 * sqrt(0.5)  -3.004261596374099e-017
    -5.282181998520569e-001 8.225260197957544e-018
    2.399888152683055e-001  -5.733259464062257e-017
    -1.209426510735405e-002 6.365366430137442e-018
    1.846294773620827e-005  -4.145515408357529e-018
    6.120179006378732e-004  -1.425483598645417e-017
    5.499197872127104e-004  3.338200835541673e-017
    3.330512243432287e-004  -3.991445470685697e-018
    1.174539197140093e-004  8.134817722064635e-018
    1.271175953190277e-003  1.061081289889354e-016
    5.904772354521236e-004  4.261935782680764e-018
    2.547449951754077e-004  -1.989653107221370e-018
    1.439335536944256e-004  1.032861020696383e-017
    8.379821248584243e-005  2.637946926618362e-018
    6.459175728329905e-005  -4.337403675070662e-018
    8.917856010668035e-004  5.921434018767376e-018];

token = reshape(token, [1, 32]);
out=zeros(1, 768);
for i = 1:768
    for j = 1:32
        out(i) = out(i) + token(j) * cos(pi*(j-1)*(2*i-1)/(768));
    end
end

out = out * sqrt(2/768) * sqrt(256) * sqrt(128)/2048;
filter_coeff=out;
tic;
total_input_f_pcm=aec('input_data/rec_mic_0_0_short.pcm', 'input_data/rec_mic_1_0_short.pcm', 'input_data/rec_spk_l_0_short.pcm', filter_coeff);
toc;
disp(['运行时间: ',num2str(toc)]);
tmp1 = total_input_f_pcm(:,1:129);
tmp2 = total_input_f_pcm(:,130:258);
[out1, out2] = compose(tmp1, tmp2, filter_coeff);
x = size(out1, 1);
y = size(out1, 2);
out1 = reshape(out1.', [x*y, 1]);
out2 = reshape(out2.', [x*y, 1]);

file_id=fopen('out_aec_0_0_debug.pcm','wb');
fwrite(file_id, out1,'int16');
fclose(file_id);

subplot(2,1,1);
plot(out1);
grid on;
subplot(2,1,2);
specgram(out1, 2048, 16000, 2048, 1024);
% % Plot the STFT result
% set(gcf, 'Position', [20 100 600 500]);
% axes('Position', [0.1 0.1 0.85 0.5]);
% n2=1:129;
% N=size(tmp1, 1);
% time=(1:N)/16000;
% freq=(n2-1)*16000/256;
% imagesc(time, freq, abs(tmp1.'));
% axis xy;ylabel('f/Hz');xlabel('Time/s');
% title('语谱图');
% m=64;
% LightYellow=[0.6 0.6 0.6];
% MidRed=[0 0 0];
% Black=[0.5 0.7 1];
% Colors=[LightYellow;MidRed;Black];
% colormap(SpecColorMap(m, Colors));
% 
% % Plot waveform
% axes('Position', [0.07 0.72 0.9 0.22]);
% plot(time, out1, 'k');
% xlim([0 max(time)]);
% xlabel('Time/s');ylabel('幅值');
% title('语音信号波形');

file_id=fopen('out_aec_1_0.pcm','wb');
fwrite(file_id, out2,'int16');
fclose(file_id);
