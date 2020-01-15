clear
clc

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

MicNum_ctrl = 6;
MicNo_ctrl = [1 2 3 4 5 6];%[1 4 6];%
MicNum_abm = 6;
MicNo_abm = [1 2 3 4 5 6];%[1 4 6];%

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
sig_mic_tao = sig_mic_tao - min(sig_mic_tao);
dist_mic_tao = dist_mic_tao - min(dist_mic_tao);
% sig_mic_tao = R_cir / c * cos(sig_angle + mic_angle);
% dist_mic_tao = R_cir / c * cos(dist_angle + mic_angle);
sig_steering = zeros(MicNum, FreLen);
dist_steering = zeros(MicNum, FreLen);
for f = 1 : FreLen
    fre = (f + FreMin-1) / FrameLen * fs;
    omiga = 2 * pi * fre;
	sig_steering(:, f) = exp(-j * omiga * sig_mic_tao);
    dist_steering(:, f) = exp(-j * omiga * dist_mic_tao);
end

%% load micArray signal
filepath = 'D:\ДњТы\gsc_yuzhou\gsc_matlab\xiaoyu_data\';


wavename = 'rec_mic_0_0.pcm';
pathInfo = dir([filepath wavename]);
SigLen = pathInfo.bytes/2;
sig_mic = zeros(MicNum, SigLen);
for m = 1 : MicNum
    wavename = ['rec_mic_' num2str(m-1) '_0.pcm'];
    pathInfo = dir([filepath wavename]);
    L = pathInfo.bytes/2;
    if L > SigLen
        L = SigLen;
    end
    fid = fopen([filepath wavename],'rb');
    sig_mic(m, 1:L) = fread(fid, [1,L], 'int16');
    fclose(fid);
end
sig_mic = sig_mic(:, 90*fs : 150*fs);

%% gsc ctrl process
tic
SigLen = size(sig_mic, 2);
FrameNum = floor((SigLen-FrameLen)/FrameInc);
gsc_delta_con = 0.0001;
gsc_delta_dyn = 0.00001;
gsc_s0_dyn = 0.00001;

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
ctrl_sig_fix = zeros(MicNum_ctrl, len128 + len32);
ctrl_dist_beam = zeros(1, len128 + len32);
ctrl_dist_fix = zeros(MicNum_ctrl, len128 + len32);
ctrl_snr_ori = zeros(FrameNum*4, len65);
ctrl_snr_pro1 = zeros(FrameNum*4, len65);
ctrl_snr_mcra = zeros(FrameNum*4, len65);
ctrl_snr_alpha = zeros(FrameNum*4, len65);
ctrl_nsr_ori = zeros(FrameNum*4, len65);
ctrl_nsr_pro1 = zeros(FrameNum*4, len65);
ctrl_nsr_mcra = zeros(FrameNum*4, len65);
ctrl_nsr_alpha = zeros(FrameNum*4, len65);
ctrl_final_alpha = zeros(FrameNum*4, len65);
ctrl_final_mcra = zeros(FrameNum*4, len65);
ctrl_abm_alpha = zeros(FrameNum*4, len65);
ctrl_aic_alpha = zeros(FrameNum*4, len65);
ctrl_abm_mcra = zeros(FrameNum*4, len65);
ctrl_aic_mcra = zeros(FrameNum*4, len65);

NoiseMcra1.first  = 1;
NoiseMcra1.len    = 65;
NoiseMcra1.w    = 1;
NoiseMcra1.L    = 60;
NoiseMcra1.alpha    = 0.9;
NoiseMcra1.frm_cnt    = 0;
NoiseMcra1.S   = zeros(1, NoiseMcra1.len);
NoiseMcra1.Smin   = zeros(1, NoiseMcra1.len);
NoiseMcra1.Stmp   = zeros(1, NoiseMcra1.len);
NoiseMcra1.Yprob   = zeros(1, NoiseMcra1.len);
NoiseMcra1.lamda_d   = zeros(1, NoiseMcra1.len);
NoiseMcra1.b   = zeros(1, 2*NoiseMcra1.w+1);    %hanning window coefficients
NoiseMcra1.b = 0.5 * (1 - cos(2 * pi * (1:2*NoiseMcra1.w+1) / (2*NoiseMcra1.w+2)));
NoiseMcra1.b = NoiseMcra1.b ./ sum(NoiseMcra1.b);

NoiseMcra2.first  = 1;
NoiseMcra2.len    = 65;
NoiseMcra2.w    = 1;
NoiseMcra2.L    = 60;
NoiseMcra2.alpha    = 0.9;
NoiseMcra2.frm_cnt    = 0;
NoiseMcra2.S   = zeros(1, NoiseMcra2.len);
NoiseMcra2.Smin   = zeros(1, NoiseMcra2.len);
NoiseMcra2.Stmp   = zeros(1, NoiseMcra2.len);
NoiseMcra2.Yprob   = zeros(1, NoiseMcra2.len);
NoiseMcra2.lamda_d   = zeros(1, NoiseMcra2.len);
NoiseMcra2.b   = zeros(1, 2*NoiseMcra2.w+1);    %hanning window coefficients
NoiseMcra2.b = 0.5 * (1 - cos(2 * pi * (1:2*NoiseMcra2.w+1) / (2*NoiseMcra2.w+2)));
NoiseMcra2.b = NoiseMcra2.b ./ sum(NoiseMcra2.b);

hwt = waitbar(0, 'GSC process');
for p = 2 : FrameNum
    waitbar(p/FrameNum, hwt);
    
    % data prepare
    ll = (p-1)*FrameInc + 1 : (p-1)*FrameInc + FrameLen;
    frame_data = sig_mic(:, ll);
    frame_data = frame_data';
    for m = 1 : MicNum
        frame_data(:,m) = frame_data(:,m) .* hamm_win;
    end
    fft_data = fft(frame_data);
    
    % signal steering and fix beam
    for m = 1 : MicNum
        fft_fix = zeros(FrameLen, 1);
        fft_fix(1 : FreLen) = fft_data(1 : FreLen, m) .* sig_steering(m, 1 : FreLen)';
        fft_fix(FrameLen/2+2:FrameLen) = conj(fft_fix(FrameLen/2:-1:2));
        fix_temp = real(ifft(fft_fix));
        fix_temp = fix_temp';
        sig_fix(m, ll) = sig_fix(m, ll) + fix_temp;
    end
    sig_beam(ll) = mean(sig_fix(MicNo_ctrl, ll));
%     sig_beam(ll) = mean(sig_fix(:, ll));
    
    % disturb steering and fix beam
    for m = 1 : MicNum
        fft_fix = zeros(FrameLen, 1);
        fft_fix(1 : FreLen) = fft_data(1 : FreLen, m) .* dist_steering(m, 1 : FreLen)';
        fft_fix(FrameLen/2+2:FrameLen) = conj(fft_fix(FrameLen/2:-1:2));
        fix_temp = real(ifft(fft_fix));
        fix_temp = fix_temp';
        dist_fix(m, ll) = dist_fix(m, ll) + fix_temp;
    end
    dist_beam(ll) = mean(dist_fix(MicNo_ctrl, ll));
%     dist_beam(ll) = mean(dist_fix(:, ll));
    
%     sig_fix(:, ll) = sig_fix_bak(:, ll);
%     dist_fix(:, ll) = dist_fix_bak(:, ll);
%     sig_beam(ll) = sig_beam_bak(ll);
%     dist_beam(ll) = dist_beam_bak(ll);
    
    for q = 1 : 4
        ll1 = (q-1)*len32 + 1 : (q-1)*len32 + len32;
        ll2 = (p-1)*FrameInc + ll1;
        
        % update control for ABM and AIC
        ctrl_sig_beam(1 : len128) = ctrl_sig_beam(len32+1 : len32+len128);
        ctrl_sig_fix(:, 1 : len128) = ctrl_sig_fix(:, len32+1 : len32+len128);
        ctrl_dist_beam(1 : len128) = ctrl_dist_beam(len32+1 : len32+len128);
        ctrl_dist_fix(:, 1 : len128) = ctrl_dist_fix(:, len32+1 : len32+len128);
        
        ctrl_sig_beam(len128+1 : len128+len32) = sig_beam(ll2);
        ctrl_sig_fix(:, len128+1 : len128+len32) = sig_fix(MicNo_ctrl, ll2);
        ctrl_dist_beam(len128+1 : len128+len32) = dist_beam(ll2);
        ctrl_dist_fix(:, len128+1 : len128+len32) = dist_fix(MicNo_ctrl, ll2);
        
        ctrl_sig_beam_fft = fft(ctrl_sig_beam(1 : len128));
        ctrl_sig_fix_fft = fft(ctrl_sig_fix(:, 1 : len128), [], 2);
        ctrl_dist_beam_fft = fft(ctrl_dist_beam(1 : len128));
        ctrl_dist_fix_fft = fft(ctrl_dist_fix(:, 1 : len128), [], 2);
        
        ctrl_sig_null_eng = zeros(1, len65);
        ctrl_dist_null_eng = zeros(1, len65);
        for m = 1 : MicNum_ctrl
            temp = ctrl_sig_beam_fft(1 : len65) - ctrl_sig_fix_fft(m, 1 : len65);
            ctrl_sig_null_eng = ctrl_sig_null_eng + abs(temp).^2 / MicNum_ctrl;
            temp = ctrl_dist_beam_fft(1 : len65) - ctrl_dist_fix_fft(m, 1 : len65);
            ctrl_dist_null_eng = ctrl_dist_null_eng + abs(temp).^2 / MicNum_ctrl;
        end
%         for m = 2 : MicNum_ctrl
%             temp = ctrl_sig_fix_fft(1, 1 : len65) - ctrl_sig_fix_fft(m, 1 : len65);
%             ctrl_sig_null_eng = ctrl_sig_null_eng + abs(temp).^2 / MicNum_ctrl;
%             temp = ctrl_dist_fix_fft(1, 1 : len65) - ctrl_dist_fix_fft(m, 1 : len65);
%             ctrl_dist_null_eng = ctrl_dist_null_eng + abs(temp).^2 / MicNum_ctrl;
%         end
        for n = 1 : len65
            if ctrl_sig_null_eng(n) < 0.001
                ctrl_sig_null_eng(n) = 0.001;
            end
            if ctrl_dist_null_eng(n) < 0.001
                ctrl_dist_null_eng(n) = 0.001;
            end
        end
        ctrl_snr_ori((p-1)*4 + q, :) = ctrl_dist_null_eng ./ ctrl_sig_null_eng;
        ctrl_nsr_ori((p-1)*4 + q, :) = ctrl_sig_null_eng ./ ctrl_dist_null_eng;
        ctrl_snr_pro1((p-1)*4 + q, :) = ctrl_snr_ori((p-1)*4 + q, :) ./ (ctrl_snr_ori((p-1)*4 + q, :) + 1.0 * ctrl_nsr_ori((p-1)*4 + q, :));
        ctrl_nsr_pro1((p-1)*4 + q, :) = ctrl_nsr_ori((p-1)*4 + q, :) ./ (1.0 * ctrl_snr_ori((p-1)*4 + q, :) + ctrl_nsr_ori((p-1)*4 + q, :));
        [ctrl_snr_mcra((p-1)*4 + q, :), NoiseMcra1] = noise_mcra(ctrl_snr_pro1((p-1)*4 + q, :), NoiseMcra1);
        [ctrl_nsr_mcra((p-1)*4 + q, :), NoiseMcra2] = noise_mcra(ctrl_nsr_pro1((p-1)*4 + q, :), NoiseMcra2);
        
        for n = 1 : len65
            if ctrl_snr_alpha((p-1)*4 + q, n) > ctrl_snr_ori((p-1)*4 + q, n)
                ctrl_snr_alpha((p-1)*4 + q, n) = 0.75 * ctrl_snr_alpha((p-1)*4 + q - 1, n) + 0.25 * ctrl_snr_pro1((p-1)*4 + q, n);
            else
                ctrl_snr_alpha((p-1)*4 + q, n) = 0.95 * ctrl_snr_alpha((p-1)*4 + q - 1, n) + 0.05 * ctrl_snr_pro1((p-1)*4 + q, n);
            end
            if ctrl_nsr_alpha((p-1)*4 + q, n) > ctrl_nsr_ori((p-1)*4 + q, n)
                ctrl_nsr_alpha((p-1)*4 + q, n) = 0.75 * ctrl_nsr_alpha((p-1)*4 + q - 1, n) + 0.25 * ctrl_nsr_pro1((p-1)*4 + q, n);
            else
                ctrl_nsr_alpha((p-1)*4 + q, n) = 0.95 * ctrl_nsr_alpha((p-1)*4 + q - 1, n) + 0.05 * ctrl_nsr_pro1((p-1)*4 + q, n);
            end
        end
        ctrl_final_alpha((p-1)*4 + q, :) = ctrl_snr_alpha((p-1)*4 + q, :) ./ (ctrl_snr_alpha((p-1)*4 + q, :) + ctrl_nsr_alpha((p-1)*4 + q, :));
%         ctrl_final_mcra((p-1)*4 + q, :) = ctrl_snr_mcra((p-1)*4 + q, :) ./ (ctrl_snr_mcra((p-1)*4 + q, :) + ctrl_nsr_mcra((p-1)*4 + q, :));
        ctrl_final_mcra((p-1)*4 + q, :) = ctrl_snr_mcra((p-1)*4 + q, :) ./ ctrl_nsr_mcra((p-1)*4 + q, :);
%         for n = 1 : len65
%             if ctrl_final_alpha((p-1)*4 + q, n) < 0.5
%                 ctrl_aic_alpha((p-1)*4 + q, n) = 1;
%             elseif ctrl_final_alpha((p-1)*4 + q, n) > 0.6
%                 ctrl_abm_alpha((p-1)*4 + q, n) = 1;
%             end
%             
%             if ctrl_final_mcra((p-1)*4 + q, n) < 0.53
%                 ctrl_aic_mcra((p-1)*4 + q, n) = 1;
%             elseif ctrl_final_mcra((p-1)*4 + q, n) > 0.63
%                 ctrl_abm_mcra((p-1)*4 + q, n) = 1;
%             end
%         end
    end
%     
    if p == 100
        p = p;
    end
end
close(hwt)
toc

%% save result
frameL = 32;
t1 = 1;
t2 = 55;
sL1 = floor(t1*fs/frameL);
sL2 = floor(t2*fs/frameL);
tt = (sL1:sL2)*frameL/fs;
ff = (1:65)/65*fs/2;
% figure,
% subplot(311), imagesc(tt, ff, ctrl_snr_pro1(sL1:sL2, :)', [0 1]), colormap(gray); title('ctrl SNR pro1');
% subplot(312), imagesc(tt, ff, ctrl_snr_alpha(sL1:sL2, :)', [0 1]), colormap(gray); title('ctrl SNR alpha');
% subplot(313), imagesc(tt, ff, ctrl_snr_mcra(sL1:sL2, :)', [0 1]), colormap(gray); title('ctrl SNR mcra');
% 
% figure,
% subplot(311), imagesc(tt, ff, ctrl_nsr_pro1(sL1:sL2, :)', [0 1]), colormap(gray); title('ctrl NSR pro1');
% subplot(312), imagesc(tt, ff, ctrl_nsr_alpha(sL1:sL2, :)', [0 1]), colormap(gray); title('ctrl NSR alpha');
% subplot(313), imagesc(tt, ff, ctrl_nsr_mcra(sL1:sL2, :)', [0 1]), colormap(gray); title('ctrl NSR mcra');
% 
% figure,
% % plot(tt, ctrl_snr_pro1(sL1:sL2, 30))
% hold on, plot(tt, ctrl_snr_alpha(sL1:sL2, 30), 'r'); grid on
% hold on, plot(tt, ctrl_snr_mcra(sL1:sL2, 30), 'c'); grid on

% figure
% subplot(311), imagesc(tt, ff, ctrl_final_alpha(sL1:sL2, :)', [0 1]), colormap(gray); title('ctrl final alpha');
% subplot(312), imagesc(tt, ff, ctrl_aic_alpha(sL1:sL2, :)'), colormap(gray); title('AIC update');
% subplot(313), imagesc(tt, ff, ctrl_abm_alpha(sL1:sL2, :)'), colormap(gray); title('ABM update');

figure
subplot(111), imagesc(tt, ff, ctrl_final_mcra(sL1:sL2, :)', [0 1]), colormap(gray); title('ctrl final mcra');
% subplot(312), imagesc(tt, ff, ctrl_aic_mcra(sL1:sL2, :)'), colormap(gray); title('AIC update');
% subplot(313), imagesc(tt, ff, ctrl_abm_mcra(sL1:sL2, :)'), colormap(gray); title('ABM update');
