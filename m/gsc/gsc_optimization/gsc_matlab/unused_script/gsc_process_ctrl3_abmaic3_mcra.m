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

MicNum_ctrl = 3;
MicNo_ctrl = [1 4 6];%[1 2 3 4 5 6];%
MicNum_abm = 3;
MicNo_abm = [1 4 6];%[1 2 3 4 5 6];%

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
filepath = 'E:\myWork\study\MicArray\DMA\xiaoyu_data\65db_3m_d1\';
% wavename = 'debug_distb_fix_0.pcm';
% pathInfo = dir([filepath wavename]);
% SigLen = pathInfo.bytes/2;
% sig_fix_bak = zeros(MicNum, SigLen);
% dist_fix_bak = zeros(MicNum, SigLen);
% for m = 1 : MicNum
%     wavename = ['debug_sig_fix_' num2str(m-1) '.pcm'];
%     pathInfo = dir([filepath wavename]);
%     L = pathInfo.bytes/2;
%     if L > SigLen
%         L = SigLen;
%     end
%     fid = fopen([filepath wavename],'rb');
%     sig_fix_bak(m, 1:L) = fread(fid, [1,L], 'int16');
%     fclose(fid);
%     
%     wavename = ['debug_distb_fix_' num2str(m-1) '.pcm'];
%     pathInfo = dir([filepath wavename]);
%     L = pathInfo.bytes/2;
%     if L > SigLen
%         L = SigLen;
%     end
%     fid = fopen([filepath wavename],'rb');
%     dist_fix_bak(m, 1:L) = fread(fid, [1,L], 'int16');
%     fclose(fid);
% end
% sig_beam_bak = mean(sig_fix_bak, 1);
% dist_beam_bak = mean(dist_fix_bak, 1);
% wavename = 'debug_ctrl_abm.pcm';
% pathInfo = dir([filepath wavename]);
% L = pathInfo.bytes/2 / 65;
% fid = fopen([filepath wavename],'rb');
% temp = fread(fid, [65, L], 'int16');
% ctrl_abm_bak = temp';
% fclose(fid);
% wavename = 'debug_ctrl_aic.pcm';
% pathInfo = dir([filepath wavename]);
% L = pathInfo.bytes/2 / 65;
% fid = fopen([filepath wavename],'rb');
% temp = fread(fid, [65, L], 'int16');
% ctrl_aic_bak = temp';
% fclose(fid);

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
ctrl_abm = zeros(FrameNum*4, len65);
ctrl_aic = zeros(FrameNum*4, len65);
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

% ABM parameter define
abm_out = zeros(MicNum_abm, SigLen);
abm_sig_beam = zeros(1, len64 + len32);
abm_sig_fix = zeros(MicNum_abm, len128);
abm_psf = zeros(MicNum_abm, len65);
abm_ht = zeros(MicNum_abm, len128);
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
aic_sig_fix = zeros(MicNum_abm, len128);
aic_psf = zeros(1, len65);
aic_hf = zeros(MicNum_abm, len65);
aic_lambda = 0.985 * (1 - 1/3/128)^(32);
aic_mu = 2 * 0.3 * (1 - aic_lambda);
aic_nu = 1 - exp(-128 / (2 * 2 * 100 * fs));
aic_norm_max = 0.02;
aic_norm_rec = zeros(FrameNum*4, 2);

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
        ctrl_final_mcra((p-1)*4 + q, :) = ctrl_snr_mcra((p-1)*4 + q, :) ./ (ctrl_snr_mcra((p-1)*4 + q, :) + ctrl_nsr_mcra((p-1)*4 + q, :));
        for n = 1 : len65
            if ctrl_final_mcra((p-1)*4 + q, n) < 0.53
                ctrl_aic((p-1)*4 + q, n) = 1;
            elseif ctrl_final_mcra((p-1)*4 + q, n) > 0.63
                ctrl_abm((p-1)*4 + q, n) = 1;
            end
        end
%         ctrl_aic((p-1)*4 + q, n) = ctrl_aic_bak((p-1)*4 + q, n);
%         ctrl_abm((p-1)*4 + q, n) = ctrl_abm_bak((p-1)*4 + q, n);
        
        % ABM process
        abm_sig_beam(1 : len64) = abm_sig_beam(len32+1 : len32+len64);
        abm_sig_beam(len64+1 : len64+len32) = mean(sig_fix(MicNo_abm, ll2));%sig_beam(ll2);%
        abm_sig_fix(:, 1 : len128-len32) = abm_sig_fix(:, len32+1 : len128);
        abm_sig_fix(:, len128-len32+1 : len128) = sig_fix(MicNo_abm, ll2);
        abm_sig_fix_fft = fft(abm_sig_fix(:, 1 : len128), [], 2);
        for m = 1 : MicNum_abm
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
        aic_sig_beam(len128-len32+1 : len128) = mean(sig_fix(MicNo_abm, ll2));%sig_beam(ll2);%
        aic_sig_fix(:, 1 : len128-len32) = aic_sig_fix(:, len32+1 : len128);
        aic_sig_fix(:, len128-len32+1 : len128) = abm_out(:, ll2);
        aic_sig_fix_fft = fft(aic_sig_fix(:, 1 : len128), [], 2);
        aic_yf = zeros(1, len128);
        aic_yf(1:len65) = sum(aic_sig_fix_fft(:, 1:len65) .* aic_hf);
        aic_yf(len65+1 : len128) = conj(aic_yf(len64:-1:2));
        temp = real(ifft(aic_yf));
        e = zeros(1, len128);
        e(len64+1 : len128) = aic_sig_beam(1 : len64) - temp(len64+1 : len128);
        ef = fft(e);
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
        for m = 1 : MicNum_abm
            aic_hf(m, :) = aic_hf(m, :) + aic_muf .* ef(1:len65) .* conj(aic_sig_fix_fft(m, 1:len65));
            aic_norm = aic_norm + sum(abs(aic_hf(m, :)).^2);
            aic_hf(m, :) = aic_hf(m, :) - aic_nuf .* aic_hf(m, :);
        end
        aic_norm = aic_norm / 128 /128;
        aic_norm_rec((p-1)*4+q, 1) = aic_norm;
        if aic_norm > aic_norm_max
            aic_norm = sqrt(aic_norm_max / aic_norm);
        else
            aic_norm = 1;
        end
        aic_norm_rec((p-1)*4+q, 2) = aic_norm;
        for m = 1 : MicNum_abm
            hf_temp = zeros(1, len128);
            hf_temp(1:len65) = aic_hf(m, 1:len65);
            hf_temp(len65+1 : len128) = conj(hf_temp(len64:-1:2));
            ht_temp = real(ifft(hf_temp));
            ht_temp(len64+1 : len128) = 0;
            ht_temp = ht_temp * aic_norm;
            hf_temp = fft(ht_temp);
            aic_hf(m, :) = hf_temp(1:len65);
        end
    end
    
    if p == 100
        p = p;
    end
end
close(hwt)
toc

%% save result
% MicNum = 1;
% for m = 1 : MicNum
%     wavename = ['rec_mic_sig_fix_' num2str(m-1) '.wav'];
%     filetemp = [filepath wavename];
%     audiowrite(filetemp, sig_fix(m,:)/32768, fs);
% end
% wavename = 'rec_mic_sig_beam.wav';
% filetemp = [filepath wavename];
% audiowrite(filetemp, sig_beam/32768, fs);
% wavename = 'rec_mic_dist_beam.wav';
% filetemp = [filepath wavename];
% audiowrite(filetemp, dist_beam/32768, fs);
% for m = 1 : MicNum_abm
%     wavename = ['rec_mic_abm_' num2str(m-1) '.wav'];
%     filetemp = [filepath wavename];
%     audiowrite(filetemp, abm_out(m,:)/32768, fs);
% end
wavename = 'rec_mic_aic.wav';
filetemp = [filepath wavename];
audiowrite(filetemp, aic_out/32768, fs);

frameL = 32;
t1 = 1;
t2 = 55;
sL1 = floor(t1*fs/frameL);
sL2 = floor(t2*fs/frameL);
tt = (sL1:sL2)*frameL/fs;
ff = (1:65)/65*fs/2;
% figure,
% subplot(211), imagesc(tt, ff, ctrl_snr_disp(sL1:sL2, :)', [0 1]), colormap(gray); title('ctrl SNRp');
% subplot(212), imagesc(tt, ff, ctrl_nsr_disp(sL1:sL2, :)', [0 1]), colormap(gray); title('ctrl NSRp');

figure, plot(tt, aic_norm_rec(sL1:sL2)), title('aic norm'); grid on

figure,
subplot(211), imagesc(tt, ff, ctrl_aic(sL1:sL2, :)'), colormap(gray); title('AIC update');
subplot(212), imagesc(tt, ff, ctrl_abm(sL1:sL2, :)'), colormap(gray); title('ABM update');

% figure,
% subplot(211), imagesc(tt, ff, ctrl_aic_bak(sL1:sL2, :)'), colormap(gray); title('AIC update');
% subplot(212), imagesc(tt, ff, ctrl_abm_bak(sL1:sL2, :)'), colormap(gray); title('ABM update');