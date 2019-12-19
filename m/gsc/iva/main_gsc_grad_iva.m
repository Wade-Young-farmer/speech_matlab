clc
clear
close all
addpath('./filter_bank');
fs = 16000;
frame_len = 8;
frame_size = fs * frame_len / 1000;
over_sample_ratio = 2;
subband_num = frame_size * over_sample_ratio / 2 + 1;
delay = 3;
h_size = frame_size * over_sample_ratio * 3;
h = filter_bank_win(h_size, delay);
filepath='../../aec/input_data/';
waveName1 = 'rec_mic_0_0_short.pcm';
waveName2 = 'rec_mic_1_0_short.pcm';
fid = fopen([filepath waveName1],'rb');
pathInfo = dir([filepath waveName1]);
L = pathInfo.bytes/2;
x1_iva = fread(fid,[L,1],'int16');
fclose(fid);
fid = fopen([filepath waveName2],'rb');
x2_iva = fread(fid,[L,1],'int16');
fclose(fid);
% x1_iva = audioread('D:\data\cube\data2\rec_mic_1_0.wav');
% x2_iva = audioread('D:\data\cube\data2\rec_mic_2_0.wav');
l2= length(x1_iva);
l1= 1;
NF = 32768;
x1_iva= x1_iva(l1:l2);
x2_iva= x2_iva(l1:l2);
nfft = 256;
win = 2*hanning(nfft,'periodic')/nfft;
nol = fix(1*2*nfft/4);
frame_num = ceil(size(x1_iva, 1) / frame_size);
S1 = filter_bank_analyze(x1_iva, frame_size, over_sample_ratio, h);%conj(stft(x1_iva, nfft, win, nol)');%
S2 = filter_bank_analyze(x2_iva, frame_size, over_sample_ratio, h);%conj(stft(x2_iva, nfft, win, nol)');%

d_mic = 0.058;
X(1,:,:) = S1.';%conj(stft(x1_iva, nfft, win, nol)');%
X(2,:,:) = S2.';%conj(stft(x2_iva, nfft, win, nol)');%
x(1,:)=x1_iva;
x(2,:)=x2_iva;
theta = 79;

[S, W,M] = gc_rtiva(X,2,theta,d_mic);
theta = 1;
MZ = 0:theta:180;
MZ = MZ/180*pi;
l = length(MZ);
nfreq = 129;
L1 = zeros(nfreq,l);
L2 = zeros(nfreq,l);
MicNum = 2;
c= 342;
fs = 16000;
for zz = 1:1:length(MZ);
    H = zeros(MicNum,nfreq);
    r=[d_mic,0];
    for f = 1 : nfreq
        fre = (f-1) / ((nfreq-1)*2) * fs;
        omiga = 2 * pi * fre;
        d0 = zeros(MicNum,1);
        for i = 1:MicNum
            d0(i) = exp(-1j * omiga*(d_mic-r(i))/c *cos(MZ(zz)));
        end
        H(:,f) = d0;
        L1(f,zz) = abs(W(1,:,f)*d0);%/(norm(W(1,:,f))*norm(d0)+1e-6);
        L2(f,zz) = abs(W(2,:,f)*d0);%/(norm(W(2,:,f))*norm(d0)+1e-6);
    end
    
end

A1 = L1./(L1+L2+1e-6);
A2 = L2./(L1+L2+1e-6);
figure
subplot(211)
imagesc(L1);
subplot(212)
imagesc(L2);
figure
plot(M(1,:),'r'),hold on;
plot(M(2,:),'b'),hold off;

X1 = squeeze(S(1,:,:)).';
X2 = squeeze(S(2,:,:)).';
psd_full = abs(S1).^2+abs(S2).^2;
psd_full  = psd_full*0.5;
psd_X1 = abs(X1).^2;
psd_X2 = abs(X2).^2;
G1 = psd_X1./(psd_full+1e-6);
G2 = psd_X2./(psd_full+1e-6);
G1 = min(G1,1);
G2 = min(G2,1);
g1_null = sqrt(1-G1).*S1;
g2_null = sqrt(1-G2).*S1;
figure
imagesc(G2);
y(1,:) = filter_bank_synthesis(X1, frame_size, over_sample_ratio, h);%istft(X1, l2, win, nol)';%
y(2,:) = filter_bank_synthesis(X2, frame_size, over_sample_ratio, h);%istft(X2, l2, win, nol)';%
y_null(1,:) = filter_bank_synthesis(g1_null, frame_size, over_sample_ratio, h);%istft(g1_null, l2, win, nol)';%
y_null(2,:) = filter_bank_synthesis(g2_null, frame_size, over_sample_ratio, h);%istft(g2_null, l2, win, nol)';%
audiowrite([filepath 'bssi_13.wav'],y(1,:)/NF*sqrt(128),16000);
audiowrite([filepath 'bssi_23.wav'],y(2,:)/NF*sqrt(128),16000);
audiowrite([filepath 'bss_1_null.wav'],y_null(1,:)/NF,16000);
audiowrite([filepath 'bss_3_null.wav'],y_null(2,:)/NF,16000);
% audiowrite('D:\data\cube\data2\bss_gc1-160-5.wav',y(1,:)/32768,16000);
% audiowrite('D:\data\cube\data2\bss_gc2-160-5.wav',y(2,:)/32768,16000);
