clear
% clc

ll = 128;
fs = 16000;
refectFlag = 0;

filepath = '/Volumes/root/baidu/personal-code/mic_array_align/icode_neon_diff_new/songruitao-3m-4F/';
waveName1 = 'opt_outsig_asr_neon.pcm';
waveName2 = 'opt_outsig_asr_icode.pcm';
pathInfo = dir([filepath waveName1]);
L = pathInfo.bytes/2;
fid = fopen([filepath waveName1],'rb');
sigMic1 = fread(fid,[L,1],'int16');
fclose(fid);
fid = fopen([filepath waveName2],'rb');
sigMic2 = fread(fid,[L,1],'int16');
fclose(fid);

figure
subplot(311); plot(sigMic1); grid on
subplot(312); plot(sigMic2); grid on
subplot(313); plot(sigMic1-sigMic2); grid on