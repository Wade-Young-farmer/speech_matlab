close all
clc
clear

% moments = [47.00811006267665, 53.05853704568263, 61.989149017868954, 68.16328923800685, 68.95018946214208, 74.24663327843685, 74.24663327843685, 80.53549708102877, 80.53549708102877, 86.94643588631027, 88.31370314267308, 93.42855459955203, 93.7614739251477, 99.63296021292591, 99.63296021292591, 105.52537042255975, 106.52412839934676, 111.51791828328183, 112.03242996829333, 117.3288737845881, 118.54413091604155, 123.44712462026871, 124.50641338352767, 129.71206102011453, 131.07400371573317, 136.55203989144377, 136.55203989144377, 141.90901449421048, 142.87250581150647, 149.47036153694793, 150.34805794079108, 155.2813170382542, 157.52095613771598, 163.7556271443258, 176.64466753941454, 180.39757630067484, 187.17702438553212, 193.13930685301824, 194.72739170366913, 199.47905844171643, 200.5383472049754, 206.22824113333778, 207.43885686277656, 213.52220090320657, 214.97258974385028, 221.38885310987595, 223.32583827697803, 228.8341398459246, 238.31848885027318, 243.00962480184856, 245.37032547425423, 252.66428524412305, 254.69206659093305, 261.4412492825544, 264.8179958150178, 271.5066477201672, 274.61666109016176, 280.7302705238277, 284.483179285088, 290.6573195052259, 291.942750755613, 297.9352986163351, 300.1446723225609, 305.5319123185636, 306.7727934412384, 311.4639293928137, 527.851643118506, 533.7590499176224, 536.3858011386621, 542.4086143926202, 547.3418734900832, 553.3041559575694, 553.576544496693, 559.4480307844713, 561.0304605905634, 566.8414160918696, 567.9612356416005, 572.6826369864118, 574.4077644008621, 579.3107581050892, 580.4305776548201, 584.9144975878546, 586.3672364631811, 591.6334148862398, 593.600665446578, 599.9563980261316, 604.5613791923305, 611.1785209261867, 612.9944445203449, 618.1395613704598, 619.1486157938202, 624.4044977703055, 625.0201020815985, 630.0138919655335, 631.1942423017363, 635.9459090397836, 637.2170555556943, 643.3306649893602, 644.1074140976335, 649.131469334462, 650.7052697700951, 656.3648982597755];
moments = [0.00000, 30.0000];
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
filepath = './test_c/';
wavename = '0.pcm';
pathInfo = dir([filepath wavename]);
SigLen = pathInfo.bytes/2;
sig_mic = zeros(MicNum, SigLen);
for m = 1 : MicNum
    wavename = [num2str(m-1) '.pcm'];
    pathInfo = dir([filepath wavename]);
    L = pathInfo.bytes/2;
    if L > SigLen
        L = SigLen;
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
gsc_delta_con = 0.0001;
gsc_delta_dyn = 0.00001;
gsc_s0_dyn = 0.00001;

%{
%% generate signal and disturb steering
R_cir = 0.03;
R_tar = 2;
sig_angle = -100 / 180 * pi;
dist_angle = -40 / 180 * pi;
mic_angle = -(0:MicNum-1)'/MicNum * 2 * pi;
mic_coord = zeros(3, MicNum);
mic_coord(1, :) = R_cir * cos(mic_angle');
mic_coord(2, :) = R_cir * sin(mic_angle');
mic_coord(1, :) = [0.03 0.015 0.015 -0.015 -0.03 -0.015];
mic_coord(2, :) = [0 -0.02598 0.02598 0.02598 0 -0.02598];
sig_coord = zeros(3, 1);
sig_coord(1) = R_tar * cos(sig_angle);
sig_coord(2) = R_tar * sin(sig_angle);
dist_coord = zeros(3, 1);
dist_coord(1) = R_tar * cos(dist_angle);
dist_coord(2) = R_tar * sin(dist_angle);
sig_mic_tao = zeros(MicNum, 1);
dist_mic_tao = zeros(MicNum, 1);
for m = 1 : MicNum
    sig_mic_tao(m) = sqrt(sum((sig_coord - mic_coord(:, m)).^2)) / c;
    dist_mic_tao(m) = sqrt(sum((dist_coord - mic_coord(:, m)).^2)) / c;
end

%% tao是距离，对齐到离声源最近的麦克风
sig_mic_tao = sig_mic_tao - min(sig_mic_tao);
dist_mic_tao = dist_mic_tao - min(dist_mic_tao);
sig_mic_tao = sig_mic_tao * fs;
dist_mic_tao = dist_mic_tao * fs;

sig_steering = zeros(MicNum, FreLen);
dist_steering = zeros(MicNum, FreLen);
for f = 1 : FreLen
    % FreMin from 1 to 129, not 0 to 128
    fre = (f + FreMin-1) / FrameLen;
    omiga = 2 * pi * fre;
    
	sig_steering(:, f) = exp(-j * omiga * sig_mic_tao);
    dist_steering(:, f) = exp(-j * omiga * dist_mic_tao);
end

% fix beam parameter define
sig_fix = zeros(MicNum, SigLen);
sig_beam = zeros(1, SigLen);
dist_fix = zeros(MicNum, SigLen);
dist_beam = zeros(1, SigLen);

% control parameter define
len128 = 128;
len65 = 65;
len64 = 64;
len32 = 32;
ctrl_sig_beam = zeros(1, len128 + len32);
ctrl_sig_fix = zeros(MicNum, len128 + len32);
ctrl_dist_beam = zeros(1, len128 + len32);
ctrl_dist_fix = zeros(MicNum, len128 + len32);

%% 一帧会有4个ctrl的值， 每个也是65个频点
ctrl_abm = zeros(FrameNum*4, len65);
ctrl_aic = zeros(FrameNum*4, len65);
ctrl_snr_disp = zeros(FrameNum*4, len65);
ctrl_nsr_disp = zeros(FrameNum*4, len65);
ctrl_snr_alpha = zeros(1, len65);
ctrl_nsr_alpha = zeros(1, len65);

% ABM parameter define
abm_out = zeros(MicNum, SigLen);
abm_sig_beam = zeros(1, len64 + len32);
abm_sig_fix = zeros(MicNum, len128);
abm_psf = zeros(MicNum, len65);
abm_ht = zeros(MicNum, len128);
abm_ht(:, 33) = 1;
temp = fft(abm_ht, [], 2);
abm_hf = temp(:, 1:len65);

abm_lambda = 0.99 * (1 - 1/3/128)^(32);
abm_mu = 2 * 0.5 * (1 - abm_lambda);
abm_nu = 1 - exp(-128 / (2 * 2 * 100 * fs));

abm_upper_bound = zeros(1, len65) + 0.001;
abm_lower_bound = zeros(1, len65) - 0.001;
abm_upper_bound(33) = 1.3;
abm_upper_bound([32, 34]) = 0.6;
abm_upper_bound([31, 35]) = 0.15;

% AIC parameter define
aic_out = zeros(1, SigLen);
aic_sig_beam = zeros(1, len128);
aic_sig_fix = zeros(MicNum, len128);
aic_psf = zeros(1, len65);
aic_hf = zeros(MicNum, len65);
aic_lambda = 0.985 * (1 - 1/3/128)^(32);
aic_mu = 2 * 0.3 * (1 - aic_lambda);
aic_nu = 1 - exp(-128 / (2 * 2 * 100 * fs));
aic_norm_max = 0.01;
%}
hwt = waitbar(0, 'GSC process');

% ATF estimate
% atf estimate parameter define
Frame_len_atf = 1024;
Frame_inc_atf = Frame_len_atf / 2;
FreLen_atf = Frame_len_atf / 2 + 1;
hamm_win_atf = hamming(Frame_len_atf, 'periodic');

FrameNum_atf = floor((SigLen-Frame_len_atf)/Frame_inc_atf);

moments_frame = floor((moments * fs - Frame_len_atf)/Frame_inc_atf); 
moment_frame_size = size(moments_frame, 2);
moment_frame_index = 1;

moment_length = 0;
moment_length_index = 0;
atf_moments_plot_1 = zeros(FreLen_atf, moment_frame_size/2);
atf_moments_plot_2 = zeros(FreLen_atf, moment_frame_size/2);
atf_moments_plot_3 = zeros(FreLen_atf, moment_frame_size/2);
atf_moments_plot_4 = zeros(FreLen_atf, moment_frame_size/2);
atf_moments_plot_5 = zeros(FreLen_atf, moment_frame_size/2);
atf_moments_plot_6 = zeros(FreLen_atf, moment_frame_size/2);


voice1 = enframe(sig_mic(1, :), hamm_win_atf, Frame_inc_atf);
voice2 = enframe(sig_mic(2, :), hamm_win_atf, Frame_inc_atf);
voice3 = enframe(sig_mic(3, :), hamm_win_atf, Frame_inc_atf);
voice4 = enframe(sig_mic(4, :), hamm_win_atf, Frame_inc_atf);
voice5 = enframe(sig_mic(5, :), hamm_win_atf, Frame_inc_atf);
voice6 = enframe(sig_mic(6, :), hamm_win_atf, Frame_inc_atf);

atf_length = 20;
h_mat_plot_1 = zeros(FreLen_atf, FrameNum_atf);
h_mat_plot_2 = zeros(FreLen_atf, FrameNum_atf);
h_mat_plot_3 = zeros(FreLen_atf, FrameNum_atf);
h_mat_plot_4 = zeros(FreLen_atf, FrameNum_atf);
h_mat_plot_5 = zeros(FreLen_atf, FrameNum_atf);
h_mat_plot_6 = zeros(FreLen_atf, FrameNum_atf);

sig_mic_reconstruct = zeros(SigLen, MicNum);
sig_mic_reconstruct = sig_mic';
%只能存储某一刻的ATF值
h_mat = ones(FreLen_atf, MicNum);

G = zeros(FreLen_atf, MicNum - 1);
G_auxiliary = zeros(FreLen_atf, MicNum - 1);
P_est = zeros(1, FreLen_atf);
%{
for p = atf_length : FrameNum_atf
    lp = (p - 1) * Frame_inc_atf + 1:(p - 1) * Frame_inc_atf + Frame_len_atf;
    
    % Get the estimate of ATF Matrix
    atf_frame_length = p - atf_length + 1:p;
    fai1 = zeros(FreLen_atf, MicNum);
    fai3 = zeros(FreLen_atf, MicNum);
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
    for j = 1: MicNum
        h_mat(:, j) = (fai1(:, j) - fai3(:, 1) .* fai3(:, j)) ./ (fai1(:, 1) - fai3(:, 1) .* fai3(:, 1));
    end
    
    % Get the output of GSC
    FBF = zeros(1, FreLen_atf);
    Y = zeros(1, FreLen_atf);
    
    for i = 1:FreLen_atf
        % FBF
        W0 = h_mat(i, :) / sum(h_mat(i, :) .* conj(h_mat(10, :)));
        W0 = conj(W0);
        FBF(i) = W0 * fv(i, :)';
        % Noise
        h_mat_block = [-conj(h_mat(i, 2:MicNum)); eye(MicNum - 1)];
        U_NOISE = conj(h_mat_block)' * fv(i, :)';
        
        % output
        Y(i) = FBF(i) - G(i, :) * U_NOISE;
    end
    
    
    % Update the G Matrix
    for i = 1: FreLen_atf
        for j = 1 : MicNum
            
        end
        
        for j = 1 : MicNum -1
            G_auxiliary(i, j) = G(i, j);
            G(i, j) = temp;
        end
    end    
end

return;
%}
for p = atf_length : FrameNum_atf
    % lp = (p - 1) * Frame_inc_atf + 1:(p - 1) * Frame_inc_atf + Frame_len_atf;
    
    if p < moments_frame(moment_frame_index)
        continue;
    elseif p > moments_frame(moment_frame_index + 1)
        moment_length_index = moment_length_index + 1;
        
        atf_moments_plot_1(:, moment_length_index) = mean(h_mat_plot_1(:, p-moment_length:p-1), 2);
        atf_moments_plot_2(:, moment_length_index) = mean(h_mat_plot_2(:, p-moment_length:p-1), 2);
        atf_moments_plot_3(:, moment_length_index) = mean(h_mat_plot_3(:, p-moment_length:p-1), 2);
        atf_moments_plot_4(:, moment_length_index) = mean(h_mat_plot_4(:, p-moment_length:p-1), 2);
        atf_moments_plot_5(:, moment_length_index) = mean(h_mat_plot_5(:, p-moment_length:p-1), 2);
        atf_moments_plot_6(:, moment_length_index) = mean(h_mat_plot_6(:, p-moment_length:p-1), 2);
        moment_length = 0;
        
        moment_frame_index = moment_frame_index + 2;
        
        if moment_frame_index > moment_frame_size
            break;
        end
        continue;
    else
        moment_length = moment_length + 1;
        %lp = (p - 1) * Frame_inc_atf + 1:(p - 1) * Frame_inc_atf + Frame_len_atf;
        atf_frame_length = p - atf_length + 1 : p;
    
        fai1 = zeros(FreLen_atf, MicNum);
        % h_mat = zeros(FreLen_atf, MicNum);
        fai3 = zeros(FreLen_atf, MicNum);
    
        out_signal = zeros(Frame_len_atf, MicNum);
    
        for i = atf_frame_length
            fv1 = fft(voice1(i, :));
            %return;
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
        
        h_mat_plot_1(:, p) = h_mat(:, 1);
        %return;
        h_mat_plot_2(:, p) = h_mat(:, 2);
        h_mat_plot_3(:, p) = h_mat(:, 3);
        h_mat_plot_4(:, p) = h_mat(:, 4);
        h_mat_plot_5(:, p) = h_mat(:, 5);
        h_mat_plot_6(:, p) = h_mat(:, 6);
        %{
        for j = 1: MicNum
            out_signal(1:FreLen_atf, j) = h_mat(:, j) .* fv(:, 1);
            out_signal(Frame_len_atf/2 + 2:Frame_len_atf, j) = conj(out_signal(Frame_len_atf/2: -1: 2, j));
            out_signal_temp = real(ifft(out_signal(1:Frame_len_atf, j)));
            sig_mic_reconstruct(lp, j) = sig_mic_reconstruct(lp, j) + out_signal_temp;
        end
        %}
    end   
end

%return;

moment_length_index = 0;
moment_frame_index = 1;
for p = atf_length : FrameNum_atf
    % lp = (p - 1) * Frame_inc_atf + 1:(p - 1) * Frame_inc_atf + Frame_len_atf;
    
    if p < moments_frame(moment_frame_index)
        continue;
    elseif p > moments_frame(moment_frame_index + 1)
        moment_length_index = moment_length_index + 1;
        moment_frame_index = moment_frame_index + 2;
        
        if moment_frame_index > moment_frame_size
            break;
        end
        continue;
    else
        lp = (p - 1) * Frame_inc_atf + 1:(p - 1) * Frame_inc_atf + Frame_len_atf;
        % atf_frame_length = p - atf_length + 1 : p;

        out_signal = zeros(Frame_len_atf, MicNum);
        fv1 = fft(voice1(p, :));
        fv1 = fv1';
        fv1 = fv1(1:FreLen_atf, :);

        out_signal(1:FreLen_atf, 1) = atf_moments_plot_1(:, moment_length_index + 1) .* fv1;
        out_signal(Frame_len_atf/2 + 2:Frame_len_atf, 1) = conj(out_signal(Frame_len_atf/2: -1: 2, 1));
        out_signal_temp = real(ifft(out_signal(1:Frame_len_atf, 1)));
        sig_mic_reconstruct(lp, 1) = sig_mic_reconstruct(lp, 1) + out_signal_temp;
        
        out_signal(1:FreLen_atf, 2) = atf_moments_plot_2(:, moment_length_index + 1) .* fv1;
        out_signal(Frame_len_atf/2 + 2:Frame_len_atf, 2) = conj(out_signal(Frame_len_atf/2: -1: 2, 2));
        out_signal_temp = real(ifft(out_signal(1:Frame_len_atf, 2)));
        sig_mic_reconstruct(lp, 2) = sig_mic_reconstruct(lp, 2) + out_signal_temp;
            
        out_signal(1:FreLen_atf, 3) = atf_moments_plot_3(:, moment_length_index + 1) .* fv1;
        out_signal(Frame_len_atf/2 + 2:Frame_len_atf, 3) = conj(out_signal(Frame_len_atf/2: -1: 2, 3));
        out_signal_temp = real(ifft(out_signal(1:Frame_len_atf, 3)));
        sig_mic_reconstruct(lp, 3) = sig_mic_reconstruct(lp, 3) + out_signal_temp;
            
        out_signal(1:FreLen_atf, 4) = atf_moments_plot_4(:, moment_length_index + 1) .* fv1;
        out_signal(Frame_len_atf/2 + 2:Frame_len_atf, 4) = conj(out_signal(Frame_len_atf/2: -1: 2, 4));
        out_signal_temp = real(ifft(out_signal(1:Frame_len_atf, 4)));
        sig_mic_reconstruct(lp, 4) = sig_mic_reconstruct(lp, 4) + out_signal_temp;
            
        out_signal(1:FreLen_atf, 5) = atf_moments_plot_5(:, moment_length_index + 1) .* fv1;
        out_signal(Frame_len_atf/2 + 2:Frame_len_atf, 5) = conj(out_signal(Frame_len_atf/2: -1: 2, 5));
        out_signal_temp = real(ifft(out_signal(1:Frame_len_atf, 5)));
        sig_mic_reconstruct(lp, 5) = sig_mic_reconstruct(lp, 5) + out_signal_temp;
            
        out_signal(1:FreLen_atf, 6) = atf_moments_plot_6(:, moment_length_index + 1) .* fv1;
        out_signal(Frame_len_atf/2 + 2:Frame_len_atf, 6) = conj(out_signal(Frame_len_atf/2: -1: 2, 6));
        out_signal_temp = real(ifft(out_signal(1:Frame_len_atf, 6)));
        sig_mic_reconstruct(lp, 6) = sig_mic_reconstruct(lp, 6) + out_signal_temp;    
    end   
end



return;

%% save result for atf
result_folder = './test_c/';
for m = 1 : MicNum
    wavname_ori = ['rec_mic_orii_' num2str(m) '.wav'];
    wavname_construct = ['rec_mic_constructt_' num2str(m) '.wav'];
    filetemp = [result_folder wavname_ori];
    audiowrite(filetemp, sig_mic(m, :)'/32768, fs);
    filetemp = [result_folder wavname_construct];
    audiowrite(filetemp, sig_mic_reconstruct(:, m)/32768, fs);
end

%return;

%***************************************
for p = 1 : FrameNum
    waitbar(p/FrameNum, hwt);
    
    % data prepare
    % 一帧是256个点
    ll = (p-1)*FrameInc + 1 : (p-1)*FrameInc + FrameLen;
    frame_data = sig_mic(:, ll);
    frame_data = frame_data';
    
    %% 分帧，分帧通过相乘，不是通过卷积
    for m = 1 : MicNum
        frame_data(:,m) = frame_data(:,m) .* hamm_win;
    end
    %% 256点的FFT
    fft_data = fft(frame_data);
    
    % signal steering and fix beam
    for m = 1 : MicNum
        fft_fix = zeros(FrameLen, 1);
        %对称性，2对应256， 128对应130， 选取1到129个点，乘以相移因子
        fft_fix(1 : FreLen) = fft_data(1 : FreLen, m) .* sig_steering(m, 1 : FreLen)';
        fft_fix(FrameLen/2+2:FrameLen) = conj(fft_fix(FrameLen/2:-1:2));
        fix_temp = real(ifft(fft_fix));
        fix_temp = fix_temp';
      
        sig_fix(m, ll) = sig_fix(m, ll) + fix_temp;
        %% 做了相移，再做了逆fft，得到的信号和原始的相差很大，没有简单的时移差，不是简单的时延
        %sig_fix(m, ll) = fix_temp;
        %return;
    end
    % fixbeam ??
    sig_beam(ll) = mean(sig_fix(:, ll));
    
    % disturb steering and fix beam
    for m = 1 : MicNum
        fft_fix = zeros(FrameLen, 1);
        fft_fix(1 : FreLen) = fft_data(1 : FreLen, m) .* dist_steering(m, 1 : FreLen)';
        fft_fix(FrameLen/2+2:FrameLen) = conj(fft_fix(FrameLen/2:-1:2));
        fix_temp = real(ifft(fft_fix));
        fix_temp = fix_temp';
        dist_fix(m, ll) = dist_fix(m, ll) + fix_temp;
        %dist_fix(m, ll) = fix_temp;
    end
    dist_beam(ll) = mean(dist_fix(:, ll));
    
%     sig_fix(:, ll) = sig_fix_bak(:, ll);
%     dist_fix(:, ll) = dist_fix_bak(:, ll);
%     sig_beam(ll) = sig_beam_bak(ll);
%     dist_beam(ll) = dist_beam_bak(ll);


    %在一帧256个点的基础上，再做处理256个点的前128个点，相当于前面的fixbeam是平滑了两帧的数据
    for q = 1 : 4
        ll1 = (q-1)*len32 + 1 : (q-1)*len32 + len32;
        ll2 = (p-1)*FrameInc + ll1;
        
        % update control for ABM and AIC
        %%移位128点
        ctrl_sig_beam(1 : len128) = ctrl_sig_beam(len32+1 : len32+len128);
        ctrl_sig_fix(:, 1 : len128) = ctrl_sig_fix(:, len32+1 : len32+len128);
        ctrl_dist_beam(1 : len128) = ctrl_dist_beam(len32+1 : len32+len128);
        ctrl_dist_fix(:, 1 : len128) = ctrl_dist_fix(:, len32+1 : len32+len128);
        
        %补入32点
        ctrl_sig_beam(len128+1 : len128+len32) = sig_beam(ll2);
        ctrl_sig_fix(:, len128+1 : len128+len32) = sig_fix(:, ll2);
        ctrl_dist_beam(len128+1 : len128+len32) = dist_beam(ll2);
        ctrl_dist_fix(:, len128+1 : len128+len32) = dist_fix(:, ll2);
        
        %做128点的fft, 相当于平滑
        ctrl_sig_beam_fft = fft(ctrl_sig_beam(1 : len128));
        ctrl_sig_fix_fft = fft(ctrl_sig_fix(:, 1 : len128), [], 2);
        ctrl_dist_beam_fft = fft(ctrl_dist_beam(1 : len128));
        ctrl_dist_fix_fft = fft(ctrl_dist_fix(:, 1 : len128), [], 2);
        
        ctrl_sig_null_eng = zeros(1, len65);
        ctrl_dist_null_eng = zeros(1, len65);
        
        %% 直接用两个null beam相比
        for m = 1 : MicNum
            temp = ctrl_sig_beam_fft(1 : len65) - ctrl_sig_fix_fft(m, 1 : len65);
            ctrl_sig_null_eng = ctrl_sig_null_eng + abs(temp).^2 / MicNum;
            temp = ctrl_dist_beam_fft(1 : len65) - ctrl_dist_fix_fft(m, 1 : len65);
            ctrl_dist_null_eng = ctrl_dist_null_eng + abs(temp).^2 / MicNum;
        end
        
        for n = 1 : len65
            if ctrl_sig_null_eng(n) < 0.001
                ctrl_sig_null_eng(n) = 0.001;
            end
            if ctrl_dist_null_eng(n) < 0.001
                ctrl_dist_null_eng(n) = 0.001;
            end
        end
        
        %% 直接用两个null beam相比
        ctrl_snr = ctrl_dist_null_eng ./ ctrl_sig_null_eng;
        ctrl_nsr = ctrl_sig_null_eng ./ ctrl_dist_null_eng;
        
        %return;
        
        
        %% 继续平滑
        for n = 1 : len65
            ctrl_snr_temp = ctrl_snr(n) / (ctrl_snr(n) + 1.5 * ctrl_nsr(n)); % ??1
            ctrl_nsr_temp = ctrl_nsr(n) / (1.5 * ctrl_snr(n) + ctrl_nsr(n)); % ??1
            if ctrl_snr_alpha(n) > ctrl_snr(n)
                ctrl_snr_alpha(n) = 0.75 * ctrl_snr_alpha(n) + 0.25 * ctrl_snr_temp;
            else
                ctrl_snr_alpha(n) = 0.95 * ctrl_snr_alpha(n) + 0.05 * ctrl_snr_temp;
            end
            if ctrl_nsr_alpha(n) > ctrl_nsr(n)
                ctrl_nsr_alpha(n) = 0.75 * ctrl_nsr_alpha(n) + 0.25 * ctrl_nsr_temp;
            else
                ctrl_nsr_alpha(n) = 0.95 * ctrl_nsr_alpha(n) + 0.05 * ctrl_nsr_temp;
            end
            
            %% 继续平滑
            snr_temp = ctrl_snr_alpha(n) / (ctrl_snr_alpha(n) + ctrl_nsr_alpha(n));
            
            %% 给一帧的四段控制矩阵赋值，每一个频点赋值
            if snr_temp < 0.5
                ctrl_aic((p-1)*4 + q, n) = 1;
            elseif snr_temp > 0.6
                ctrl_abm((p-1)*4 + q, n) = 1;
            end
        end
        
        %% 画图用的中间值？
        ctrl_snr_disp((p-1)*4 + q, :) = ctrl_snr_alpha;
        ctrl_nsr_disp((p-1)*4 + q, :) = ctrl_nsr_alpha;
%         ctrl_aic((p-1)*4 + q, n) = ctrl_aic_bak((p-1)*4 + q, n);
%         ctrl_abm((p-1)*4 + q, n) = ctrl_abm_bak((p-1)*4 + q, n);



        % ABM process
        % ????32?, buffer???96
        abm_sig_beam(1 : len64) = abm_sig_beam(len32+1 : len32+len64);
        abm_sig_beam(len64+1 : len64+len32) = sig_beam(ll2);
        % ????
        abm_sig_fix(:, 1 : len128-len32) = abm_sig_fix(:, len32+1 : len128);
        % ?????
        abm_sig_fix(:, len128-len32+1 : len128) = sig_fix(:, ll2);
        
        %% ???????? 128??fft
        %{
        abm_sig_beam_fft = fft(abm_sig_beam(1:len64));
        for m = 1 : MicNum
            abm_estimator_f = zeros(1, len128);
            abm_estimator_f(1 : len65) = abm_sig_beam_fft(1:len65) .* abm_hf(m, :);
            abm_estimator_f(len65 + 1 : len128) = conj(abm_estimator_f(len64:-1:2));
            temp = real(ifft(abm_estimator_f));
            e = zeros(1, len128);
            e(len64 + 1 : len128) = abm_sig_fix(m, 1:len64) - temp(len64+1:len128);
            ef = fft(e);
        end
        %}
        abm_sig_fix_fft = fft(abm_sig_fix(:, 1 : len128), [], 2);
        for m = 1 : MicNum
            abm_yf = zeros(1, len128);
            abm_yf(1:len65) = abm_sig_fix_fft(m, 1:len65) .* abm_hf(m, :);
            abm_yf(len65+1 : len128) = conj(abm_yf(len64:-1:2));
            temp = real(ifft(abm_yf));
            e = zeros(1, len128);
            e(len64+1 : len128) = abm_sig_beam(1 : len64) - temp(len64+1 : len128);
            ef = fft(e);
            abm_out(m, ll2) = e(len128-len32+1 : len128);
            
            abm_sig_fix_fft_eng = abs(abm_sig_fix_fft(m, 1:len65)).^2;
            abm_psf(m, :) = abm_lambda * abm_psf(m, :) + (1 - abm_lambda) * abm_sig_fix_fft_eng;
            abm_psf_temp = abm_psf(m, :) .* (abm_psf(m, :) > gsc_delta_con) + gsc_delta_con * (abm_psf(m, :) <= gsc_delta_con);
            abm_muf = abm_mu * ctrl_abm((p-1)*4 + q, 1:len65) ./ abm_psf_temp;
%             if (p-1)*4 + q > 2
%                 abm_nuf = abm_nu * ctrl_aic((p-1)*4 + q - 1, 1:len65);
%             else
%                 abm_nuf = zeros(1, len65);
%             end
            abm_nuf = abm_nu * ctrl_aic((p-1)*4 + q, 1:len65);
            abm_hf(m, :) = abm_hf(m, :) + abm_muf .* ef(1:len65) .* conj(abm_sig_fix_fft(m, 1:len65));
            abm_hf(m, :) = abm_hf(m, :) - abm_nuf .* abm_hf(m, :);
            
            
            hf_temp = zeros(1, len128);
            hf_temp(1:len65) = abm_hf(m, 1:len65);
            hf_temp(len65+1 : len128) = conj(hf_temp(len64:-1:2));
            ht_temp = real(ifft(hf_temp));
            ht_temp(len64+1 : len128) = 0;
            for n = 1 : 30
                ht_temp(n) = min([ht_temp(n), abm_upper_bound(n)]);
                ht_temp(n) = max([ht_temp(n), abm_lower_bound(n)]);
            end
            for n = 36 : 64
                ht_temp(n) = min([ht_temp(n), abm_upper_bound(n)]);
                ht_temp(n) = max([ht_temp(n), abm_lower_bound(n)]);
            end
            hf_temp = fft(ht_temp);
            abm_hf(m, :) = hf_temp(1:len65);
        end
        aa = 1;
        
        % AIC process
        aic_sig_beam(1 : len128-len32) = aic_sig_beam(len32+1 : len128);
        % ?abm?????fixedbeam??????????
        aic_sig_beam(len128-len32+1 : len128) = sig_beam(ll2);
        % ???????32?abm???
        aic_sig_fix(:, 1 : len128-len32) = aic_sig_fix(:, len32+1 : len128);
        aic_sig_fix(:, len128-len32+1 : len128) = abm_out(:, ll2);
        % ???????
        aic_sig_fix_fft = fft(aic_sig_fix(:, 1 : len128), [], 2);
        aic_yf = zeros(1, len128);
        aic_yf(1:len65) = sum(aic_sig_fix_fft(:, 1:len65) .* aic_hf);
        aic_yf(len65+1 : len128) = conj(aic_yf(len64:-1:2));
        temp = real(ifft(aic_yf));
        e = zeros(1, len128);
        % ????????????fft???64???????????(?????
        e(len64+1 : len128) = aic_sig_beam(1 : len64) - temp(len64+1 : len128);
        ef = fft(e);
        % ??? 64???32???????
        aic_out(ll2) = e(len128-len32+1 : len128);
        
        aic_sig_fix_fft_eng = abs(aic_sig_fix_fft(:, 1:len65)).^2;
        aic_sig_fix_fft_eng = sum(aic_sig_fix_fft_eng);
        aic_psf = aic_lambda * aic_psf + (1 - aic_lambda) * aic_sig_fix_fft_eng;
        aic_psf_temp = zeros(1, len65);
        for n = 1 : len65
            aic_psf_temp(n) = aic_psf(n) + gsc_delta_dyn * exp(-aic_psf(n) / gsc_s0_dyn);
            if aic_psf_temp(n) < 10e-6
                aic_psf_temp(n) = 10e-6;
            end
        end
%         if (p-1)*4 + q > 2
%             aic_muf = aic_mu * ctrl_aic((p-1)*4 + q - 1, 1:len65) ./ aic_psf_temp;
%         else
%         	aic_muf = zeros(1, len65);
%         end
        aic_muf = aic_mu * ctrl_aic((p-1)*4 + q, 1:len65) ./ aic_psf_temp;
        aic_nuf = aic_nu * ctrl_abm((p-1)*4 + q, 1:len65);
        aic_norm = 0;
        for m = 1 : MicNum
            % ???????????????
            aic_hf(m, :) = aic_hf(m, :) + aic_muf .* ef(1:len65) .* conj(aic_sig_fix_fft(m, 1:len65));
            aic_norm = aic_norm + sum(abs(aic_hf(m, :)).^2);
            aic_hf(m, :) = aic_hf(m, :) - aic_nuf .* aic_hf(m, :);
        end
        
        aic_norm = aic_norm / 128 /128;
        if aic_norm > aic_norm_max
            aic_norm = sqrt(aic_norm_max / aic_norm);
        else
            aic_norm = 1;
        end
        for m = 1 : MicNum
            hf_temp = zeros(1, len128);
            hf_temp(1:len65) = aic_hf(m, 1:len65);
            hf_temp(len65+1 : len128) = conj(hf_temp(len64:-1:2));
            ht_temp = real(ifft(hf_temp));
            % ?????????64????????L = 64, ??????????????????
            ht_temp(len64+1 : len128) = 0;
            ht_temp = ht_temp * aic_norm;
            hf_temp = fft(ht_temp);
            aic_hf(m, :) = hf_temp(1:len65);
        end
    end
    
    %{
    if p == 100
        p = p;
    end
    %}
end
close(hwt)
toc

%% save result
for m = 1 : MicNum
    wavename = ['rec_mic_sig_fix_' num2str(m-1) '.wav'];
    filetemp = [filepath wavename];
    % ???mic????????
    audiowrite(filetemp, sig_fix(m,:)/32768, fs);
end
wavename = 'rec_mic_sig_beam.wav';
filetemp = [filepath wavename];
audiowrite(filetemp, sig_beam/32768, fs);
wavename = 'rec_mic_dist_beam.wav';
filetemp = [filepath wavename];
audiowrite(filetemp, dist_beam/32768, fs);
for m = 1 : MicNum
    wavename = ['rec_mic_abm_' num2str(m-1) '.wav'];
    filetemp = [filepath wavename];
    audiowrite(filetemp, abm_out(m,:)/32768, fs);
end
wavename = 'rec_mic_aic.wav';
filetemp = [filepath wavename];
audiowrite(filetemp, aic_out/32768, fs);

frameL = 32;
t1 = 10;
t2 = 55;
sL1 = floor(t1*fs/frameL);
sL2 = floor(t2*fs/frameL);
tt = (sL1:sL2)*frameL/fs;
ff = (1:65)/65*fs/2;
% figure,
% subplot(211), imagesc(tt, ff, ctrl_snr_disp(sL1:sL2, :)', [0 1]), colormap(gray); title('ctrl SNRp');
% subplot(212), imagesc(tt, ff, ctrl_nsr_disp(sL1:sL2, :)', [0 1]), colormap(gray); title('ctrl NSRp');

figure,
subplot(211), imagesc(tt, ff, ctrl_aic(sL1:sL2, :)'), colormap(gray); title('AIC update');
subplot(212), imagesc(tt, ff, ctrl_abm(sL1:sL2, :)'), colormap(gray); title('ABM update');

% figure,
% subplot(211), imagesc(tt, ff, ctrl_aic_bak(sL1:sL2, :)'), colormap(gray); title('AIC update');
% subplot(212), imagesc(tt, ff, ctrl_abm_bak(sL1:sL2, :)'), colormap(gray); title('ABM update');