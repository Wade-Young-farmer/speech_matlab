close all
clear
clc

load ./3mic_xiaoyu/abm.txt
load ./3mic_xiaoyu/aic.txt
load ./3mic_xiaoyu/beta.txt
load ./3mic_xiaoyu/c_beta.txt

load ./3mic_xiaoyu/f_m_pf_nsr.txt
load ./3mic_xiaoyu/f_m_pf_nsr_noi.txt
load ./3mic_xiaoyu/f_m_pf_snr.txt
load ./3mic_xiaoyu/f_m_pf_snr_noi.txt

fs = 16000;
frameL = 32;
%% 时间是一分40， 到2分25
t1 = 81;
t2 = 95;
sL1 = floor(t1*fs/frameL);
sL2 = floor(t2*fs/frameL);
tt = (sL1:sL2)*frameL/fs;
ff = (1:65)/65*fs/2;
% figure,
% subplot(211), imagesc(tt, ff, ctrl_snr_disp(sL1:sL2, :)', [0 1]), colormap(gray); title('ctrl SNRp');
% subplot(212), imagesc(tt, ff, ctrl_nsr_disp(sL1:sL2, :)', [0 1]), colormap(gray); title('ctrl NSRp');

figure,
% 白色是1， 黑色是0
subplot(211), imagesc(tt, ff, beta(sL1:sL2, 1:65)'), colormap(gray); title('SNR update');
subplot(212), imagesc(tt, ff, c_beta(sL1:sL2, 1:65)'), colormap(gray); title('NSR update');

figure,
% 白色是1， 黑色是0
subplot(211), imagesc(tt, ff, abm(sL1:sL2, 1:65)'), colormap(gray); title('abm update');
subplot(212), imagesc(tt, ff, aic(sL1:sL2, 1:65)'), colormap(gray); title('aic update');

%{
figure,
% 白色是1， 黑色是0
subplot(211), imagesc(tt, ff, f_m_pf_nsr(sL1:sL2, 1:65)'), colormap(gray); title('f_m_pf_nsr update');
subplot(212), imagesc(tt, ff, f_m_pf_snr(sL1:sL2, 1:65)'), colormap(gray); title('f_m_pf_snr update');

figure,
% 白色是1， 黑色是0
subplot(211), imagesc(tt, ff, f_m_pf_nsr_noi(sL1:sL2, 1:65)'), colormap(gray); title('f_m_pf_nsr_noi update');
subplot(212), imagesc(tt, ff, f_m_pf_snr_noi(sL1:sL2, 1:65)'), colormap(gray); title('f_m_pf_snr_noi update');
%}
