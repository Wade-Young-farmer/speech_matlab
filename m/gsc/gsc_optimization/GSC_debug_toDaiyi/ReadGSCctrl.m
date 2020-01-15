clear
clc

filepath = '/Users/daiyi02/Downloads/GSC_debug_toDaiyi/bd_spil_gsc_debug/';
N = 65;
fileTemp = [filepath 'debug_gsc_ctrl_NSR.bin'];
pathInfo = dir(fileTemp);
Len = pathInfo.bytes/4/N;
L1 = floor(Len);
fid = fopen(fileTemp,'rb');
temp = fread(fid, [N, L1], 'float32');
fclose(fid);
size(temp)
ctrl_NSR = temp(:, 1:4:L1);
ctrl_NSRn = temp(:, 2:4:L1);
ctrl_NSRp = temp(:, 3:4:L1);
ctrl_ABM = temp(:, 4:4:L1);
L1 = floor(L1/4);
ctrl_NSRpsd = ctrl_NSRp(:, 1:L1) ./ ctrl_NSRn(:, 1:L1);


fileTemp = [filepath 'debug_gsc_ctrl_SNR.bin'];
pathInfo = dir(fileTemp);
Len = pathInfo.bytes/4/N;
L2 = floor(Len);
fid = fopen(fileTemp,'rb');
temp = fread(fid, [N, L2], 'float32');
fclose(fid);
size(temp)
ctrl_SNR = temp(:, 1:4:L2);
ctrl_SNRn = temp(:, 2:4:L2);
ctrl_SNRp = temp(:, 3:4:L2);
ctrl_AIC = temp(:, 4:4:L2);
L2 = floor(L2/4);
ctrl_SNRpsd = ctrl_SNRp(:, 1:L2) ./ ctrl_SNRn(:, 1:L2);


ctrl_SNRsnr = ctrl_SNRp./ctrl_NSRp;

fs = 16000;
frameL = 32;
% t1 = 1;
% t2 = min(L1, L2) * (frameL / fs);
t1 = 55;
t2 = 75;
% t1 = 60;
% t2 = 75 ;
sL1 = floor(t1*fs/frameL);
sL2 = floor(t2*fs/frameL);
tt = (sL1:sL2)*frameL/fs;
ff = (1:N)/N*fs/2;

% figure,
% subplot(211), imagesc(tt, ff, ctrl_SNR(:,sL1:sL2), [0 10]), colormap(gray); title('ctrl SNR');
% subplot(212), imagesc(tt, ff, ctrl_NSR(:,sL1:sL2), [0 10]), colormap(gray); title('ctrl NSR');

% figure,
% subplot(211), imagesc(tt, ff, ctrl_SNRp(:,sL1:sL2), [0 10]), colormap(gray); title('ctrl SNRp');
% subplot(212), imagesc(tt, ff, ctrl_NSRp(:,sL1:sL2), [0 10]), colormap(gray); title('ctrl NSRp');

figure,
subplot(211), imagesc(tt, ff, ctrl_SNRsnr(:,sL1:sL2)), colormap(jet); title('ctrl SNRsnr');
subplot(212), imagesc(tt, ff, ctrl_SNRp(:,sL1:sL2)), colormap(jet); title('ctrl SNRp');

% figure,
% subplot(211), imagesc(tt, ff, ctrl_SNRn(:,sL1:sL2), [0 10]), colormap(gray); title('ctrl SNRn');
% subplot(212), imagesc(tt, ff, ctrl_NSRn(:,sL1:sL2), [0 10]), colormap(gray); title('ctrl NSRn');

% figure,
% subplot(211), imagesc(tt, ff, ctrl_SNRpsd(:,sL1:sL2), [0 30]), colormap(gray); title('ctrl SNRpsd');
% subplot(212), imagesc(tt, ff, ctrl_NSRpsd(:,sL1:sL2), [0 30]), colormap(gray); title('ctrl NSRpsd');

% figure,
% subplot(211), imagesc(tt, ff, ctrl_AIC(:,sL1:sL2)), colormap(gray); title('ctrl AIC');
% subplot(212), imagesc(tt, ff, ctrl_ABM(:,sL1:sL2)), colormap(gray); title('ctrl ABM');

% figure
% subplot(211), plot(tt, sum(ctrl_SNRn(:,sL1:sL2))), title('ctrl AIC'); grid on
% subplot(212), plot(tt, sum(ctrl_NSRn(:,sL1:sL2))), title('ctrl ABM'); grid on

% figure
% subplot(211), plot(tt, ctrl_SNRn(1,sL1:sL2), 'b'), title('ctrl AIC'); grid on
% hold on, plot(tt, ctrl_SNRn(2,sL1:sL2), 'r'), title('ctrl AIC'); grid on
% subplot(212), plot(tt, ctrl_NSRn(1,sL1:sL2)), title('ctrl ABM'); grid on
% hold on, plot(tt, ctrl_SNRn(2,sL1:sL2), 'r'), title('ctrl ABM'); grid on

% figure
% subplot(211), plot(tt, ctrl_SNRn(3,sL1:sL2)), title('ctrl AIC'); grid on;axis([min(tt) max(tt) 0 10])
% subplot(212), plot(tt, ctrl_NSRn(3,sL1:sL2)), title('ctrl ABM'); grid on;axis([min(tt) max(tt) 0 10])

% figure
% plot(tt, ctrl_SNR(10, sL1:sL2), 'b'); grid on
% hold on, plot(tt, ctrl_SNRp(5, sL1:sL2), 'r');grid on
% hold on, plot(tt, ctrl_SNRn(5, sL1:sL2), 'g');grid on
% % hold on, plot(tt, ctrl_SNRpsd(5, sL1:sL2)/10, 'm');grid on
% axis([min(tt) max(tt) 0 10])


% temp1 = ctrl_SNRp(:,sL1:sL2) < 2.0;
% temp2 = ctrl_NSRp(:,sL1:sL2) < 1.2;
% temp = temp1 .* temp2;
% figure,
% subplot(211), imagesc(tt, ff, temp1 - temp), colormap(gray); title('SNR');
% subplot(212), imagesc(tt, ff, temp2 - temp), colormap(gray); title('NSR');

% temp1 = ctrl_SNRsnr(:,sL1:sL2) < 0.5;           % AIC update
% temp2 = ctrl_SNRsnr(:,sL1:sL2) > 4.0;           % ABM update
% temp = temp1 .* temp2;
% figure,
% subplot(211), imagesc(tt, ff, temp1 - temp), colormap(gray); title('AIC update');
% subplot(212), imagesc(tt, ff, temp2 - temp), colormap(gray); title('ABM update');

figure,
subplot(211), imagesc(tt, ff, ctrl_AIC(:,sL1:sL2)), colormap(gray); title('AIC update');
subplot(212), imagesc(tt, ff, ctrl_ABM(:,sL1:sL2)), colormap(gray); title('ABM update');

% tdsip = 167.8;
% Ldisp = floor(tdsip*fs/frameL);
% figure,
% subplot(211), plot(ctrl_SNRp(:,Ldisp)); grid on
% subplot(212), plot(ctrl_NSRp(:,Ldisp)); grid on
% 
% temp1 = reshape(ctrl_SNRn(:,sL1:sL2), 1, N*(sL2-sL1+1));
% temp2 = reshape(ctrl_NSRn(:,sL1:sL2), 1, N*(sL2-sL1+1));
% [SNRa, SNRb] = hist(temp1, 0:0.01:100);
% [NSRa, NSRb] = hist(temp2, 0:0.01:100);
% figure, 
% subplot(211);plot(SNRb, SNRa); title('SNR'); grid on;axis([min(SNRb) max(SNRb) 0 100])
% subplot(212);plot(NSRb, NSRa); title('NSR'); grid on;axis([min(NSRb) max(NSRb) 0 100])

