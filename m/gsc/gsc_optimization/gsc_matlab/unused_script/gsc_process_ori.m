close all
clear

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

%% generate signal and disturb steering
R_cir = 0.03;
R_tar = 2;
sig_angle = -90 / 180 * pi;
dist_angle = -50 / 180 * pi;
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
filepath = '/Volumes/daiyi/20180705_xiaoyu_asr/185_0330so_nooneshot_wz65db/65db_1m_d1/';
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
%% 时间是一分半到两分半
sig_mic = sig_mic(:, 60*fs + 1: 120*fs);

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
ctrl_sig_fix = zeros(MicNum, len128 + len32);
ctrl_dist_beam = zeros(1, len128 + len32);
ctrl_dist_fix = zeros(MicNum, len128 + len32);
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

hwt = waitbar(0, 'GSC process');
for p = 1 : FrameNum
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
    end
    dist_beam(ll) = mean(dist_fix(:, ll));
    
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
        ctrl_sig_fix(:, len128+1 : len128+len32) = sig_fix(:, ll2);
        ctrl_dist_beam(len128+1 : len128+len32) = dist_beam(ll2);
        ctrl_dist_fix(:, len128+1 : len128+len32) = dist_fix(:, ll2);
        
        ctrl_sig_beam_fft = fft(ctrl_sig_beam(1 : len128));
        ctrl_sig_fix_fft = fft(ctrl_sig_fix(:, 1 : len128), [], 2);
        ctrl_dist_beam_fft = fft(ctrl_dist_beam(1 : len128));
        ctrl_dist_fix_fft = fft(ctrl_dist_fix(:, 1 : len128), [], 2);
        
        ctrl_sig_null_eng = zeros(1, len65);
        ctrl_dist_null_eng = zeros(1, len65);
        %% ?????????
        for m = 1 : MicNum
            %% ?????????ctrl_sig_fix_fft????????????????????????
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
        ctrl_snr = ctrl_dist_null_eng ./ ctrl_sig_null_eng;
        ctrl_nsr = ctrl_sig_null_eng ./ ctrl_dist_null_eng;
        for n = 1 : len65
            ctrl_snr_temp = ctrl_snr(n) / (ctrl_snr(n) + 1.5 * ctrl_nsr(n));
            ctrl_nsr_temp = ctrl_nsr(n) / (1.5 * ctrl_snr(n) + ctrl_nsr(n));
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
            snr_temp = ctrl_snr_alpha(n) / (ctrl_snr_alpha(n) + ctrl_nsr_alpha(n));
            if snr_temp < 0.5
                ctrl_aic((p-1)*4 + q, n) = 1;
            elseif snr_temp > 0.6
                ctrl_abm((p-1)*4 + q, n) = 1;
            end
        end
        
        % ???? 
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
%{
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
%}
frameL = 32;
%% 时间是一分40， 到2分25
t1 = 21;
t2 = 35;
sL1 = floor(t1*fs/frameL);
sL2 = floor(t2*fs/frameL);
tt = (sL1:sL2)*frameL/fs;
ff = (1:65)/65*fs/2;
% figure,
% subplot(211), imagesc(tt, ff, ctrl_snr_disp(sL1:sL2, :)', [0 1]), colormap(gray); title('ctrl SNRp');
% subplot(212), imagesc(tt, ff, ctrl_nsr_disp(sL1:sL2, :)', [0 1]), colormap(gray); title('ctrl NSRp');

figure,
% 白色是1， 黑色是0
subplot(212), imagesc(tt, ff, ctrl_abm(sL1:sL2, :)'), colormap(gray); title('ABM update');
subplot(211), imagesc(tt, ff, ctrl_aic(sL1:sL2, :)'), colormap(gray); title('AIC update');


% figure,
% subplot(211), imagesc(tt, ff, ctrl_aic_bak(sL1:sL2, :)'), colormap(gray); title('AIC update');
% subplot(212), imagesc(tt, ff, ctrl_abm_bak(sL1:sL2, :)'), colormap(gray); title('ABM update');