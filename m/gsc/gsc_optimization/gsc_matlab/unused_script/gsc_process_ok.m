close all
clc
clear

moments = [47.00811006267665, 53.05853704568263, 61.989149017868954, 68.16328923800685, 68.95018946214208, 74.24663327843685, 74.24663327843685, 80.53549708102877, 80.53549708102877, 86.94643588631027, 88.31370314267308, 93.42855459955203, 93.7614739251477, 99.63296021292591, 99.63296021292591, 105.52537042255975, 106.52412839934676, 111.51791828328183, 112.03242996829333, 117.3288737845881, 118.54413091604155, 123.44712462026871, 124.50641338352767, 129.71206102011453, 131.07400371573317, 136.55203989144377, 136.55203989144377, 141.90901449421048, 142.87250581150647, 149.47036153694793, 150.34805794079108, 155.2813170382542, 157.52095613771598, 163.7556271443258, 176.64466753941454, 180.39757630067484, 187.17702438553212, 193.13930685301824, 194.72739170366913, 199.47905844171643, 200.5383472049754, 206.22824113333778, 207.43885686277656, 213.52220090320657, 214.97258974385028, 221.38885310987595, 223.32583827697803, 228.8341398459246, 238.31848885027318, 243.00962480184856, 245.37032547425423, 252.66428524412305, 254.69206659093305, 261.4412492825544, 264.8179958150178, 271.5066477201672, 274.61666109016176, 280.7302705238277, 284.483179285088, 290.6573195052259, 291.942750755613, 297.9352986163351, 300.1446723225609, 305.5319123185636, 306.7727934412384, 311.4639293928137, 527.851643118506, 533.7590499176224, 536.3858011386621, 542.4086143926202, 547.3418734900832, 553.3041559575694, 553.576544496693, 559.4480307844713, 561.0304605905634, 566.8414160918696, 567.9612356416005, 572.6826369864118, 574.4077644008621, 579.3107581050892, 580.4305776548201, 584.9144975878546, 586.3672364631811, 591.6334148862398, 593.600665446578, 599.9563980261316, 604.5613791923305, 611.1785209261867, 612.9944445203449, 618.1395613704598, 619.1486157938202, 624.4044977703055, 625.0201020815985, 630.0138919655335, 631.1942423017363, 635.9459090397836, 637.2170555556943, 643.3306649893602, 644.1074140976335, 649.131469334462, 650.7052697700951, 656.3648982597755];
%% parameter define
fs = 16000;
c = 340;
MicNum = 6;
FrameLen = fs * 0.016;
FrameInc = FrameLen/2;
FreMin = floor(0 / fs * FrameLen) + 1;
FreMax = floor(fs/2 / fs * FrameLen) + 1;
FreLen = FreMax - FreMin + 1;
hamm_win = hamming(FrameLen);
fre_disp = (FreMin : FreMax) / FrameLen * fs;

%% load micArray signal
filepath = './aj_185_1m_d2/';
wavename = 'rec_mic_0_0.pcm';
pathInfo = dir([filepath wavename]);
SigLen = pathInfo.bytes/2;
sig_mic = zeros(MicNum, SigLen);
for m = 1 : MicNum
    wavename = ['rec_mic_' num2str(m-1) '_0.pcm'];
    pathInfo = dir([filepath wavename]);
    L = pathInfo.bytes/2;
    if L > SigLen
        L = SigLen;gs
    end
    fid = fopen([filepath wavename],'rb');
    sig_mic(m, 1:L) = fread(fid, [1,L], 'int16');
    fclose(fid);
end

%sig_mic = sig_mic(:, 60*fs + 1 : 120*fs);
%% gsc ctrl process
tic
SigLen = size(sig_mic, 2);
FrameNum = floor((SigLen-FrameLen)/FrameInc);
hwt = waitbar(0, 'GSC process');

% ATF estimate
% atf estimate parameter define
Frame_len_atf = 1024;
Frame_inc_atf = Frame_len_atf / 2;
FreLen_atf = Frame_len_atf / 2 + 1;
hamm_win_atf = hamming(Frame_len_atf, 'periodic');

FrameNum_atf = floor((SigLen-Frame_len_atf)/Frame_inc_atf);

voice1 = enframe(sig_mic(1, :), hamm_win_atf, Frame_inc_atf);
voice2 = enframe(sig_mic(2, :), hamm_win_atf, Frame_inc_atf);
voice3 = enframe(sig_mic(3, :), hamm_win_atf, Frame_inc_atf);
voice4 = enframe(sig_mic(4, :), hamm_win_atf, Frame_inc_atf);
voice5 = enframe(sig_mic(5, :), hamm_win_atf, Frame_inc_atf);
voice6 = enframe(sig_mic(6, :), hamm_win_atf, Frame_inc_atf);

atf_length = 20;
h_mat_plot_2 = zeros(FreLen_atf, FrameNum_atf);
h_mat_plot_3 = zeros(FreLen_atf, FrameNum_atf);

sig_mic_reconstruct = zeros(SigLen, MicNum);
%只能存储某一刻的ATF值
h_mat = ones(FreLen_atf, MicNum);

G = zeros(FreLen_atf, MicNum - 1);
G_auxiliary = zeros(FreLen_atf, MicNum - 1);
P_est = zeros(1, FreLen_atf);

for p = atf_length : FrameNum_atf
    lp = (p - 1) * Frame_inc_atf + 1:(p - 1) * Frame_inc_atf + Frame_len_atf;
    atf_frame_length = p - atf_length + 1 : p;
    
    fai1 = zeros(FreLen_atf, MicNum);
     % h_mat = zeros(FreLen_atf, MicNum);
    fai3 = zeros(FreLen_atf, MicNum);
    
    out_signal = zeros(Frame_len_atf, MicNum);
    
    for i = atf_frame_length
       fv1 = fft(voice1(i, :));
       fv2 = fft(voice2(i, :));
       fv3 = fft(voice3(i, :));
       fv4 = fft(voice4(i, :));
       fv5 = fft(voice5(i, :));
       fv6 = fft(voice6(i, :));
       
       fv = [fv1;fv2;fv3;fv4;fv5;fv6]';
       fv = fv(1:FreLen_atf,:); 
       for j = 1 : MicNum
           fai3(:,j) = fai3(:, j) + fv(:, 1) .* conj(fv(:, j));
           fai1(:,j) = fai1(:, j) + fv(:, 1) .* conj(fv(:, 1)) .* fv(:, 1) .* conj(fv(:, j));
       end
    end
    %% 计算atf 矩阵
    
    %% 一次最小二乘解
    for j = 1: MicNum
        h_mat(:, j) = (fai1(:, j) - fai3(:, 1) .* fai3(:, j)) ./ (fai1(:, 1) - fai3(:, 1) .* fai3(:, 1));
    end
    
    h_mat_plot_2(:, p) = h_mat(:, 2);
    h_mat_plot_3(:, p) = h_mat(:, 3);
    
    for j = 1: MicNum
        out_signal(1:FreLen_atf, j) = h_mat(:, j) .* fv(:, 1);
        out_signal(Frame_len_atf/2 + 2:Frame_len_atf, j) = conj(out_signal(Frame_len_atf/2: -1: 2, j));
        out_signal_temp = real(ifft(out_signal(1:Frame_len_atf, j)));
        sig_mic_reconstruct(lp, j) = sig_mic_reconstruct(lp, j) + out_signal_temp;
    end  
end

%return;

%% save result for atf
result_folder = './aj_185_1m_d2/';
for m = 1 : MicNum
    wavname_ori = ['rec_mic_ori_' num2str(m) '.wav'];
    wavname_construct = ['rec_mic_construct_' num2str(m) '.wav'];
    filetemp = [result_folder wavname_ori];
    audiowrite(filetemp, sig_mic(m, :)'/32768, fs);
    filetemp = [result_folder wavname_construct];
    audiowrite(filetemp, sig_mic_reconstruct(:, m)/32768, fs);
end

toc
close(hwt)

return;