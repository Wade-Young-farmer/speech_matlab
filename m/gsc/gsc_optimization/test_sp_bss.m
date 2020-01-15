clear
clc
close all

fs = 16000;
c = 343;

 MicNum = 2;
 r_cir = 0.01;
 delta = 2 * pi * r_cir / MicNum;
 micAngle = (0:MicNum-1)'/MicNum * 2 * pi;
 tao01 = 2 * pi * r_cir / MicNum / c;
 PMic = zeros(3, MicNum);
 PMic(1, :) = r_cir * cos(micAngle');
 PMic(2, :) = r_cir * sin(micAngle');
nfft = 512;
freMin = floor(0 / fs * nfft)+1;
freMax = floor(fs/2 / fs * nfft);
% BeamPattern define
BeamTheta = (1:360)/180*pi;
BeamThetaNum = length(BeamTheta);
BeamPattern = zeros(nfft/2+1,BeamThetaNum);
BeamPatternc = zeros(nfft/2+1,BeamThetaNum);
BeamPatternDB = zeros(nfft/2+1,BeamThetaNum);
load top_filter_24k.txt
load end_filter_24k.txt
top_filter_16k = resample(top_filter_24k,16000,24000);
end_filter_16k = resample(end_filter_24k,16000,24000);
F1 = fft(top_filter_16k,nfft);
F2 = fft(end_filter_16k,nfft);
AF=0.8;

%% DMA simulation

DMA_Steering = zeros( nfft/2 +1,MicNum);
for f = 2 : nfft/2+1
    fre = (f + freMin-1) / nfft * fs;
    omiga = 2 * pi * fre;
    GAMMA = zeros(MicNum, MicNum) + 1;
    for m1 = 1 : MicNum
        for m2 = m1+1 : MicNum
            xtemp = omiga * r_cir * abs(sin(pi*(m1-m2)/MicNum)) / c;
            GAMMA(m1, m2) = sin(xtemp) / xtemp;
            GAMMA(m2, m1) = GAMMA(m1, m2);
        end
    end
    d0 = exp(-j * omiga * r_cir' / c .* cos(micAngle'));
    d0 = d0';
    if(f < 6)
        dn_temp = (GAMMA + 0.00 * diag(zeros(1,MicNum)+1))^(-1);
    else
        dn_temp = (GAMMA + 0.000 * diag(zeros(1,MicNum)+1))^(-1);
    end
    if(f < 1)
         DMA_Steering(f,:) =  (dn_temp * d0) / ( d0' * dn_temp * d0);
    else
        %DMA_Steering(f,:) = [F1(f),F2(f)];
        DMA_Steering(f,:) = [1/4+1i/(2*omiga/c * 0.02),1/4-1i/(2*omiga/c*0.02)];
    end
    for n = 1 : BeamThetaNum
        D = exp(j * omiga * r_cir / c .* cos(BeamTheta(n) - micAngle));
        BeamPatternc(f,n) = (D'*DMA_Steering(f,:).');
        BeamPattern(f,n) = abs(D'*DMA_Steering(f,:).');
        BeamPatternDB(f,n) = 20*log10(abs(BeamPattern(f,n)));
        
    end
    
end
FrameLen = fs * 0.016;
FrameInc = FrameLen/2;
FreDisp = [250, 500, 1000, 2000, 4000, 5500, 7100, 8000];
FreDispNo = floor(FreDisp / fs * FrameLen);
Beampatt = zeros(length(FreDisp), BeamThetaNum);
figure
dispNum = length(FreDisp)/2;
for n = 1 : length(FreDisp)
    if (n==1)
        Beampatt(n,:) = mean(BeamPattern(1:FreDispNo(n), :),1);
        subplot(2,dispNum,n),polarLz(BeamTheta, 20*log10(Beampatt(n,:)),'-b',-40,0); title(['0-' num2str(FreDisp(n)) 'Hz'])
    else
        Beampatt(n,:) = mean(BeamPattern(FreDispNo(n-1):FreDispNo(n), :),1);
        subplot(2,dispNum,n),polarLz(BeamTheta, 20*log10(Beampatt(n,:)),'-b',-40,0); title([num2str(FreDisp(n-1)) '-' num2str(FreDisp(n)) 'Hz'])
    end
end
Beampatt = 20*log10(Beampatt);

% Plot Beam patterns in the same figure
Beampatt = Beampatt(:,[BeamThetaNum/2+1:BeamThetaNum,1:BeamThetaNum/2]);
figure;
plot(BeamTheta*180/pi, Beampatt(1, :), 'k-', BeamTheta*180/pi, Beampatt(2, :), 'g-',...
        BeamTheta*180/pi, Beampatt(3, :), 'r-', BeamTheta*180/pi, Beampatt(4, :), 'm-',...
        BeamTheta*180/pi, Beampatt(5, :), 'b-', BeamTheta*180/pi, Beampatt(6, :), 'y-', 'LineWidth', 1.5);
temp = min(min(Beampatt));
temp = max(floor(temp/10)*10, -60);
axis([0 360 temp 0]);
% axis([0 360 0 1]);
h = legend('OctaveBand (0-250Hz)', 'OctaveBand (250-500Hz)', 'OctaveBand (500-1000Hz)', 'OctaveBand (1000-2000Hz)', ...
    'OctaveBand (2000-4000Hz)', 'OctaveBand (4000-8000Hz)', 'Location', 'SouthEast');
set(h, 'Fontsize', 10, 'FontWeight', 'demi');
xlabel('\it Degree', 'Fontsize', 14);
ylabel('\it Beam pattern (dB)', 'Fontsize', 14);
grid on;
%%process
if 1
    file_dir =  'D:\data\车载录音_0801\houdingdeng_ganrao_ningjun_fujiahouWakjiashihouDis\';
    fs = 16000;
    file_name1 = [file_dir,'rec_mic_6_0.wav'];
    file_name2=[file_dir,'rec_mic_1_0.wav'];
    file_name3=[file_dir,'rec_mic_3_0.wav'];
    [x_mic1,fs_in]=audioread(file_name1);
    [x_mic2,fs_in]=audioread(file_name2);
    [x_mic3,fs_in]=audioread(file_name3);
    x_mic1_in = resample(x_mic1,fs,fs_in,20);
    x_mic2_in = resample(x_mic2,fs,fs_in,20);
    x_mic3_in = resample(x_mic3,fs,fs_in,20);
    audiowrite('x_fin1.wav',x_mic1_in,fs);
    % audiowrite('x_fin2.wav',x_mic2_in,fs);
    % audiowrite('x_fin3.wav',x_mic3_in,fs);
    if 1
        X1 = stft(x_mic1_in,nfft,0.25*nfft,1);
        X2 = stft(x_mic2_in,nfft,0.25*nfft,1);
        X3 = stft(x_mic3_in,nfft,0.25*nfft,1);
        [m,n]= size(X1);
        Xout1 = zeros(m,n);
        Xout2 = zeros(m,n);
        Xout3 = zeros(m,n);
        Xout1b = zeros(m,n);
        Xout2b = zeros(m,n);
        Xout3b = zeros(m,n);
        Xout4b = zeros(m,n);
        Xtmp = zeros(m,n);
        filter1 = zeros(m,n);
        filter2 = zeros(m,n);
        for k = 1:1:n
         X1f = X1(:,k).* DMA_Steering(:,1);
         X2f = X2(:,k).* DMA_Steering(:,2);
         X3f = X3(:,k).* DMA_Steering(:,2);
         w = 1:1:nfft/2+1;
         w = 2*pi*w'/nfft*fs*0.02/c;  

         Xout1(:,k) = X1f+X2f;
         Xout2(:,k) = X1f+X3f;
         X1b = X1(:,k).* DMA_Steering(:,2);
         X2b = X2(:,k).* DMA_Steering(:,1);
         X3b = X3(:,k).* DMA_Steering(:,1);

         Xout1b(:,k) = X2b + X1b;
         Xout2b(:,k) = X3b + X1b;
         
         Xout3b(:,k) = abs(abs(Xout1b(:,k)./(1-AF*exp(-1i*w)))+abs(Xout1(:,k)./(1-AF*exp(-1i*w))) - abs(X1(:,k) + X2(:,k)));
         
         tmp = (abs(Xout1(:,k))./(abs(Xout1(:,k))+abs(Xout3b(:,k))));
         Xout1(:,k) = Xout1(:,k).*tmp;
         Xout4b(:,k) = abs(abs(Xout2b(:,k)./(1-AF*exp(-1i*w)))+abs(Xout2(:,k)./(1-AF*exp(-1i*w))) - abs(X1(:,k) + X3(:,k)));
         Xout2(:,k) = Xout2(:,k).*(abs(Xout2(:,k))./(abs(Xout2(:,k))+abs(Xout4b(:,k))));
        end
        xout1 = istft(Xout1,nfft,0.25*nfft,1);
        xout2 = istft(Xout2,nfft,0.25*nfft,1);
        xout3 = istft(Xout3b,nfft,0.25*nfft,1);
        xout4 = istft(Xout4b,nfft,0.25*nfft,1);
    end

    audiowrite('x_out1_endfire2.wav',xout1*5,fs);
    audiowrite('x_out2_endfire.wav',xout2*5,fs);
    audiowrite('x_out3_endfire.wav',xout3*5,fs);
    audiowrite('x_out4_endfire.wav',xout4*5,fs);
    %%%%%%%%%%%%%%%%%%%%%%%%
    %multistage winner filter
    if 1
        Xoutfilter = Xout1;
        InterCanFilter = zeros(nfft/2+1,1);
        SigCanFilter = zeros(nfft/2+1,1);
        psd_main = abs(Xout1).^2;%%主信号支路
        psd_ref = abs(Xout2).^2;%%噪音参考
        psd_sig=zeros(nfft/2+1,n);
        psd_filter_out = zeros(nfft/2+1,n);
        winner_filter = zeros(nfft/2+1,n);
        psd_main_smooth = zeros(nfft/2+1,1);
        psd_inter_smooth = zeros(nfft/2+1,1);
        psd_sig_smooth = zeros(nfft/2+1,1);
        ctrl_snr_smooth = zeros(nfft/2+1,1);
        ctrl_nsr_smooth = zeros(nfft/2+1,1);
        ctrl_sf = [0.75*ones(nfft/2+1,1),0.95*ones(nfft/2+1,1)];
        ctrl_sig_update = zeros(nfft/2+1,1);
        ctrl_inter_update = zeros(nfft/2+1,1);
        psd_sf = 0.0;
        %%time loop
        for t = 1:1:n
        %%preprocess:calculate the filter update control
            ctrl_snr = psd_main(:,t)./(psd_ref(:,t)+1e-6);
            ctrl_nsr = psd_ref(:,t)./(psd_main(:,t)+1e-6);
            
            for k =1:1:nfft/2+1
                if ctrl_snr_smooth(k) > ctrl_snr(k)
                    ctrl_snr_smooth(k) = 0.75 * ctrl_snr_smooth(k) + 0.25 * ctrl_snr(k);
                else
                    ctrl_snr_smooth(k) = 0.95 * ctrl_snr_smooth(k) + 0.05 * ctrl_snr(k);
                end
                if ctrl_nsr_smooth(k) > ctrl_nsr(k)
                    ctrl_nsr_smooth(k) = 0.75 * ctrl_nsr_smooth(k) + 0.25 * ctrl_nsr(k);
                else
                    ctrl_nsr_smooth(k) = 0.95 * ctrl_nsr_smooth(k) + 0.05 * ctrl_nsr(k);
                end
            end
            
            ctrl_temp = ctrl_snr_smooth ./ (ctrl_nsr_smooth + ctrl_snr_smooth + 1e-6);
            ctrl_sig_update = (ctrl_temp > 0.6);
            ctrl_inter_update = (ctrl_temp < 0.8);
        %%step1 filter the input data
            psd_sig(:,t)=psd_main(:,t)-InterCanFilter.*psd_ref(:,t); 
            psd_sig(:,t) = max(psd_sig(:,t),0);
            psd_sig(:,t) = min(psd_sig(:,t),psd_main(:,t));
        %%step2 calculate the signal psd and the winner filter gain
            winner_filter = psd_sig(:,t)./(psd_main(:,t)+1e-6);
            winner_filter = max(winner_filter,0.1);
            winner_filter = min(winner_filter,1);
            winner_filter(1) =1;
            psd_filter_out(:,t) = winner_filter.*psd_sig(:,t);
        %%step3 update the filter coeff.
            %%step3.1 smmoth the input psd
            psd_main_smooth = psd_sf*psd_main_smooth+(1-psd_sf)*psd_main(:,t);
            psd_inter_smooth = psd_sf*psd_inter_smooth+(1-psd_sf)*psd_ref(:,t);
            psd_sig_smooth = psd_sf*psd_main_smooth+(1-psd_sf)*psd_filter_out(:,t);
            %%step3.2 calcultate the err.
            psd_err = psd_main_smooth - psd_inter_smooth.*InterCanFilter-psd_sig_smooth.*SigCanFilter;
            %%step3.3 use nlms alg. to update the coeff.
            InterCanFilter = InterCanFilter+0.1*(psd_err.*ctrl_inter_update ./(psd_sig_smooth.^2 + psd_inter_smooth.^2)+1e-6).*psd_ref(:,t);
            SigCanFilter = SigCanFilter+0.1*(psd_err.*ctrl_sig_update ./(psd_sig_smooth.^2 + psd_inter_smooth.^2)+1e-6).*psd_filter_out(:,t);
        end
        Xoutfilter = sqrt(psd_filter_out).*exp(1i*angle(Xout1));
        xout1f = istft(Xoutfilter,nfft,0.25*nfft,1);
        audiowrite('x_out1_endfire_filter2.wav',xout1f,fs);
    end
    %% gsc ctrl process
if 0
tic
SigLen = length(xout2);
sig_mic = zeros(MicNum, SigLen);
sig_mic(1,:) = xout1'*32768;
sig_mic(2,:) = xout2'*32768;
FrameLen = fs * 0.016*2;
FrameInc = FrameLen/2;
freMin = floor(0 / fs * nfft)+1;
freMax = floor(fs/2 / fs * nfft);
FrameNum = floor((SigLen-FrameLen)/FrameInc);
gsc_delta_con = 0.0001;
gsc_delta_dyn = 0.00001;
gsc_s0_dyn = 0.00001;

sig_fix = zeros(MicNum, SigLen);

% fix beam parameter define
sig_beam = zeros(1, SigLen);
dist_beam = zeros(1, SigLen);
%beamsteering 
sig_mic_tao = zeros(MicNum, 1);
sig_mic_tao(1,1) = -0.02/c;%2*0.01*sin(135/180*pi/2)*sin(135/2/180*pi) / c;
sig_mic_tao(2,1) = 0;%2*0.01*sin(135/180*pi/2)*sin(135/180*pi) / c;
FreMin = floor(0 / fs * FrameLen) + 1;
FreMax = floor(fs/2 / fs * FrameLen) + 1;
FreLen = FreMax - FreMin + 1;
sig_steering = zeros(MicNum, FreLen);
for f = 1 : FreLen
    fre = (f + FreMin-1) / FrameLen * fs;
    omiga = 2 * pi * fre;
	sig_steering(:, f) = exp(-j * omiga * sig_mic_tao);
end
% control parameter define
len128 = 128*2;
len64 = 64*2;
len65 = len64+1;

len32 = 32*2;
len16 = 16*2;
ctrl_sig_beam = zeros(1, len128 + len32);
ctrl_dist_beam = zeros(1, len128 + len32);
ctrl_abm = zeros(FrameNum*4, len65);
ctrl_aic = zeros(FrameNum*4, len65);
ctrl_snr_disp = zeros(FrameNum*4, len65);
ctrl_nsr_disp = zeros(FrameNum*4, len65);
ctrl_snr_alpha = zeros(1, len65);
ctrl_nsr_alpha = zeros(1, len65);

snr_out =  zeros(1, SigLen);
% ABM parameter define
abm_out = zeros(1, SigLen);
abm_sig_beam = zeros(1, len64 + len32);
abm_sig_fix = zeros(1, len128);
abm_psf = zeros(1, len65);
abm_ht = zeros(1, len128);
abm_ht(1, len32+1) = 1;
temp = fft(abm_ht);
abm_hf = temp(1, 1:len65);
abm_lambda = 0.99 * (1 - 1/3/len128)^(len32);
abm_mu = 2 * 0.5 * (1 - abm_lambda);
abm_nu = 1 - exp(-128 / (2 * 2 * 100 * fs));
abm_upper_bound = zeros(1, len65) + 0.001;
abm_lower_bound = zeros(1, len65) - 0.001;
abm_upper_bound(len32+1) = 1.3;
abm_upper_bound([len32, len32+2]) = 0.6;
abm_upper_bound([len32-1, len32+3]) = 0.15;

% AIC parameter define
syncdly=len32+len32+4;
aic_out = zeros(1, SigLen);
aic_sig_beam = zeros(1, len64+syncdly);
aic_sig_fix = zeros(1, len128);
aic_psf = zeros(1, len65);
aic_hf = zeros(1, len65);
aic_lambda = 0.985 * (1 - 1/3/len128)^(len32);
aic_mu = 2 * 0.5 * (1 - aic_lambda);
aic_nu = 1 - exp(-len128 / (2 * 2 * 100 * fs));
aic_norm_max = 1;%0.01;

hwt = waitbar(0, 'GSC process');
hamm_win = hamming(FrameLen);
%step3 calculate the control matrix to update the block matrix and
%interference matrix
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
    sig_beam(ll) = sig_mic(1, ll);
    dist_beam(ll) = sig_mic(2, ll);
    
    
    
    for q = 1 : 4
        ll1 = (q-1)*len32 + 1 : (q-1)*len32 + len32;
        ll2 = (p-1)*FrameInc + ll1;
        
        % update control for ABM and AIC
        ctrl_sig_beam(1 : len128) = ctrl_sig_beam(len32+1 : len32+len128);
        ctrl_dist_beam(1 : len128) = ctrl_dist_beam(len32+1 : len32+len128);   
        ctrl_sig_beam(len128+1 : len128+len32) = sig_beam(ll2);
        ctrl_dist_beam(len128+1 : len128+len32) = dist_beam(ll2);
        
        ctrl_sig_beam_fft = fft(ctrl_sig_beam(1 : len128));
        ctrl_dist_beam_fft = fft(ctrl_dist_beam(1 : len128));        
        ctrl_sig_null_eng = zeros(1, len65);
        ctrl_dist_null_eng = zeros(1, len65);
        ctrl_dist_null_eng = abs(ctrl_sig_beam_fft);
        ctrl_sig_null_eng = abs(ctrl_dist_beam_fft);
        ctrl_snr = ctrl_dist_null_eng ./ (ctrl_sig_null_eng+1e-6);
        ctrl_nsr = ctrl_sig_null_eng ./ (ctrl_dist_null_eng+1e-6);
        for n = 1 : len65
            ctrl_snr_temp = ctrl_snr(n);
            ctrl_nsr_temp = ctrl_nsr(n);
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
            snr_temp = ctrl_snr_alpha(n) / (ctrl_snr_alpha(n) + ctrl_nsr_alpha(n) + 1e-6);
            if snr_temp < 0.4
                ctrl_aic((p-1)*4 + q, n) = 1;
            elseif snr_temp > 0.6
                ctrl_abm((p-1)*4 + q, n) = 1;
            end
        end
        temp_ctrl = zeros(1, len128); 
        temp_ctrl(1:len65) = ctrl_sig_beam_fft(1, 1:len65).*snr_temp;
        temp_ctrl(len65+1 : len128) = conj(temp_ctrl(len64:-1:2));
        temp = real(ifft(temp_ctrl));
        snr_out(1, ll2) = temp(len128-len32+1 : len128);
        ctrl_snr_disp((p-1)*4 + q, :) = ctrl_snr_alpha;
        ctrl_nsr_disp((p-1)*4 + q, :) = ctrl_nsr_alpha;
        if 1
            % ABM process
            %buffer process 
            abm_sig_beam(1 : len64) = abm_sig_beam(len32+1 : len32+len64);
            abm_sig_beam(len64+1 : len64+len32) = dist_beam(1, ll2);%sig_beam(ll2);
            abm_sig_fix(1, 1 : len128-len32) = abm_sig_fix(1, len32+1 : len128);
            abm_sig_fix(1, len128-len32+1 : len128) = sig_beam(1,ll2);%dist_beam(1, ll2);
            abm_sig_fix_fft = fft(abm_sig_fix(1, 1 : len128));
            %filter process        
            abm_yf = zeros(1, len128); 
            abm_yf(1:len65) = abm_sig_fix_fft(1:len65) .* abm_hf(1, :);
            abm_yf(len65+1 : len128) = conj(abm_yf(len64:-1:2));
            temp = real(ifft(abm_yf));
            e = zeros(1, len128);
            e(len64+1 : len128) = abm_sig_beam(1 : len64) - temp(len64+1 : len128);
            ef = fft(e);
            abm_out(1, ll2) = e(len128-len32+1 : len128);

            abm_sig_fix_fft_eng = abs(abm_sig_fix_fft(1:len65)).^2;
            abm_psf(1, :) = abm_lambda * abm_psf(1, :) + (1 - abm_lambda) * abm_sig_fix_fft_eng;
            abm_psf_temp = abm_psf(1, :);% .* (abm_psf(1, :) > gsc_delta_con) + gsc_delta_con * (abm_psf(1, :) <= gsc_delta_con);
            abm_muf = abm_mu * ctrl_abm((p-1)*4 + q, 1:len65) ./ abm_psf_temp;
             abm_nuf = abm_nu * ctrl_aic((p-1)*4 + q, 1:len65);
            abm_hf(1, :) = abm_hf(1, :) + abm_muf .* ef(1:len65) .* conj(abm_sig_fix_fft(1:len65));
              abm_hf(1, :) = abm_hf(1, :) - abm_nuf .* abm_hf(1, :);
             hf_temp = zeros(1, len128);
             hf_temp(1:len65) = abm_hf(1, 1:len65);
             hf_temp(len65+1 : len128) = conj(hf_temp(len64:-1:2));
             ht_temp = real(ifft(hf_temp));
             ht_temp(len64+1 : len128) = 0;
             for n = 1 : len32 -2%30
                 ht_temp(n) = min([ht_temp(n), abm_upper_bound(n)]);
                 ht_temp(n) = max([ht_temp(n), abm_lower_bound(n)]);
             end
             for n = len32 + 4 : len64
                 ht_temp(n) = min([ht_temp(n), abm_upper_bound(n)]);
                 ht_temp(n) = max([ht_temp(n), abm_lower_bound(n)]);
             end
             hf_temp = fft(ht_temp);
    %           abm_hf(1, :) = hf_temp(1:len65);
             aa = 1;
              syncdly1 = syncdly;
        else
            syncdly1 = syncdly - len32;
            abm_out(1, ll2) = dist_beam(1, ll2);
        end
        
       if 1
            % AIC process
            aic_sig_beam(1 : len64+syncdly1-len32) = aic_sig_beam(len32+1 : len64+syncdly1);
            aic_sig_beam(len64+syncdly1-len32+1 : len64+syncdly1) = sig_beam(ll2);
            aic_sig_fix(1, 1 : len128-len32) = aic_sig_fix(1, len32+1 : len128);
            aic_sig_fix(1, len128-len32+1 : len128) = abm_out(1, ll2);
            aic_sig_fix_fft = fft(aic_sig_fix(1, 1 : len128));
            aic_yf = zeros(1, len128);
            aic_yf(1:len65) = (aic_sig_fix_fft(1, 1:len65) .* aic_hf);
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
                aic_psf_temp(n) = aic_psf(n);% + gsc_delta_dyn * exp(-aic_psf(n) / gsc_s0_dyn);
                if aic_psf_temp(n) < 10e-6
                    aic_psf_temp(n) = 10e-6;
                end
            end
            aic_muf = aic_mu * ctrl_aic((p-1)*4 + q, 1:len65) ./ aic_psf_temp;
            aic_nuf = aic_nu * ctrl_abm((p-1)*4 + q, 1:len65);
            aic_norm = 0;
            aic_hf(1, :) = aic_hf(1, :) + aic_muf .* ef(1:len65) .* conj(aic_sig_fix_fft(1, 1:len65));
            aic_norm = aic_norm + sum(abs(aic_hf(1, :)).^2);
            %aic_hf(1, :) = aic_hf(1, :) - aic_nuf .* aic_hf(1, :);
            aic_norm = aic_norm / len128 /len128;
            if aic_norm > aic_norm_max
                aic_norm = sqrt(aic_norm_max / aic_norm);
            else
                aic_norm = 1;
            end
            hf_temp = zeros(1, len128);
            hf_temp(1:len65) = aic_hf(1, 1:len65);
            hf_temp(len65+1 : len128) = conj(hf_temp(len64:-1:2));
            ht_temp = real(ifft(hf_temp));
            ht_temp(len64+1 : len128) = 0;
            ht_temp = ht_temp * aic_norm;
            hf_temp = fft(ht_temp);
    %          aic_hf(1, :) = hf_temp(1:len65);
       else
           %len16 update aic filter
           loopnum = len32/len16;
           sig_buf1=sig_beam(ll2);
           dist_buf1 = abm_out(1, ll2);
           aic_out1=zeros(1,len32);
           for lp = 1:1:loopnum
               ll9 = (lp-1)*len16+1:lp*len16;
                aic_sig_beam(1 : len64+syncdly1-len16) = aic_sig_beam(len16+1 : len64+syncdly1);
                aic_sig_beam(len64+syncdly1-len16+1 : len64+syncdly1) = sig_buf1(ll9);
                aic_sig_fix(1, 1 : len128-len16) = aic_sig_fix(1, len16+1 : len128);
                aic_sig_fix(1, len128-len16+1 : len128) = dist_buf1(ll9);
                aic_sig_fix_fft = fft(aic_sig_fix(1, 1 : len128));
                %filter process
                aic_yf = zeros(1, len128);
                aic_yf(1:len65) = (aic_sig_fix_fft(1:len65) .* aic_hf);
                aic_yf(len65+1 : len128) = conj(aic_yf(len64:-1:2));
                temp = real(ifft(aic_yf));
                e = zeros(1, len128);
                e(len64+1 : len128) = aic_sig_beam(1 : len64) - temp(len64+1 : len128);
                ef = fft(e);
                aic_out1(ll9) =  e(len128-len16+1 : len128);
                %aic_out(ll2) = e(len128-len16+1 : len128);
                aic_sig_fix_fft_eng = abs(aic_sig_fix_fft(:, 1:len65)).^2;
                aic_sig_fix_fft_eng = sum(aic_sig_fix_fft_eng);
                aic_psf = abm_lambda * aic_psf + (1 - abm_lambda) * aic_sig_fix_fft_eng;
                aic_psf_temp = aic_psf;
                aic_muf = aic_mu * ctrl_aic((p-1)*4 + q, 1:len65) ./ aic_psf_temp;
                aic_hf(1, :) = aic_hf(1, :) + aic_muf .* ef(1:len65) .* conj(aic_sig_fix_fft(1, 1:len65));
           end
           aic_out(ll2)=aic_out1;
  
       end
    end
    
    if p == 100
        p = p;
    end    
end
close(hwt)
wavename = 'x_out_gsc_endfire.wav';
audiowrite(wavename, aic_out/32768, fs);
wavename = 'x_out_gsc_abm_endfire.wav';
audiowrite(wavename, abm_out/32768, fs);
frameL = 32*2;
t1 = 80;
t2 = 170;
sL1 = floor(t1*fs/frameL);
sL2 = floor(t2*fs/frameL);
tt = (sL1:sL2)*frameL/fs;
ff = (1:65)/65*fs/2;
figure,
subplot(211), imagesc(tt, ff, ctrl_aic(sL1:sL2, :)'), colormap(gray); title('AIC update');
subplot(212), imagesc(tt, ff, ctrl_abm(sL1:sL2, :)'), colormap(gray); title('ABM update');
end

    
    %%%%%%%%%%%%%%%%%%%%%%%
    %bss alg.
    if 0
        %%bss param. calculate
        thate = 90;
        d = 0.02;

        tao = d/c*sin(thate/180*pi);
        Xout2_angle = (angle(Xout2.*conj(Xout1)));
        tao_D = zeros(m,n);
        C = ones(m,n);
        sigma=max(abs(d/c*sin((thate+90)*pi/180)-tao),abs(d/c*sin((thate-90)*pi/180)-tao));
        for k = 2:1:m
            Xout2_angle(k,:) = Xout2_angle(k,:)*nfft/fs/2/pi/(k-1);
            tao_D(k,:)=abs(Xout2_angle(k,:)-tao);
        end 
        C = (tao_D < sigma);
        figure
        imagesc(C);
        Xout2 = Xout2 .* (C);
        xout21 = istft(Xout2,nfft,0.25*nfft,1);
        audiowrite('x_out2_endfire_bss1.wav',xout21,fs);
    end
end