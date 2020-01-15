clear
clc
fs = 16000;

NoiseMcra1.first  = 1;
NoiseMcra1.len    = 129;
NoiseMcra1.w    = 1;
NoiseMcra1.L    = 60;
NoiseMcra1.alpha    = 0.85;
NoiseMcra1.frm_cnt    = 0;
NoiseMcra1.S   = zeros(1, NoiseMcra1.len);
NoiseMcra1.Smin   = zeros(1, NoiseMcra1.len);
NoiseMcra1.Stmp   = zeros(1, NoiseMcra1.len);
NoiseMcra1.Yprob   = zeros(1, NoiseMcra1.len);
NoiseMcra1.lamda_d   = zeros(1, NoiseMcra1.len);
NoiseMcra1.b   = zeros(1, 2*NoiseMcra1.w+1);    %hanning window coefficients
NoiseMcra1.b = 0.5 * (1 - cos(2 * pi * (1:2*NoiseMcra1.w+1) / (2*NoiseMcra1.w+2)));
NoiseMcra1.b = NoiseMcra1.b ./ sum(NoiseMcra1.b);

alafa0 = 2;
s = 1/10;  

energy_data = zeros(1,NoiseMcra1.len);

fout1 = fopen('D:\ДњТы\result_1.pcm','wb');
fout2 = fopen('D:\ДњТы\result_2.pcm','wb');

FrameLen = 256;
FrameInc = FrameLen / 2;
hamm_win = hamming(FrameLen);

noi_spec = zeros(1, NoiseMcra1.len);
lastsig = zeros(1, FrameLen/2);

fid = fopen('out_left_fix.pcm','rb');
sig_mic = fread(fid, 'int16');
fclose(fid);


SigLen = length(sig_mic);
FrameNum = floor((SigLen - FrameLen) / FrameInc);
f_sig = zeros(NoiseMcra1.len, FrameNum);
for p = 1 : FrameNum
    ls = (p - 1) *  FrameInc + 1 : (p - 1) *  FrameInc + 256;
    tmp_sig_mic(:,p) = sig_mic(ls);
end
% testsig = enframe(sig_mic,hanning(FrameLen,'periodic'),FrameInc);
%testsig = testsig';
%ftestsig = fft(testsig);
%testfft_data = abs(ftestsig(1:129,:));

fft_data = (stft(sig_mic,FrameLen,0.5 * FrameLen,1));

hwt = waitbar(0, 'spe process');
for p = 1 : FrameNum
    waitbar(p/FrameNum, hwt);
    
    energy_data = 0.8 * energy_data + 0.2 * abs(fft_data(:,p))';
    
    [noi_spec, NoiseMcra1] = mcra(abs(fft_data(:,p))',NoiseMcra1);
    
    snr = 10 * log10((energy_data .^ 2) ./ (noi_spec .^ 2));
    
    alafa = alafa0 - s * snr;
    alafa = max(alafa,1);
    alafa = min(alafa,3);
    
    yizhi(:,p) = max((energy_data - alafa .* noi_spec),0) ./ (energy_data + eps);
    yizhi(:,p) = max(yizhi(:,p),0.01);
    
%     for i = 1 : 5
%         if(p > 1)
%             yizhi(i,p) = 0.8 * yizhi(i,p - 1) + 0.2 * yizhi(i,p);
%         end
%     end

    f_sig(: ,p) = (yizhi(:,p) .* fft_data(:,p));
%     
%     track_noi(:,p - 6910) = noi_spec;
%     real_noi(:,p - 6910) = energy_data;
%     pro(:,p - 6910) = real_noi(:,p - 6910) ./ track_noi(:,p - 6910);
    
    track_noi(:,p) = noi_spec;
    real_noi(:,p) = energy_data;
    pro(:,p) = real_noi(:,p) ./ track_noi(:,p);
    

%     f_sig(FrameLen / 2 + 2 : FrameLen) = fliplr(f_sig(2 : FrameLen / 2));
%     tsig = real(ifft(f_sig));
%     outsig = tsig(1:FrameLen / 2) + lastsig;
%     final = int16(outsig);
% %     fprintf(fout1,'%d\r\n',final);
%     fwrite(fout1,final,'int16');
%     fwrite(fout2,int16(noi_spec),'int16');
%     lastsig = tsig(FrameLen / 2 + 1 : FrameLen);
end
    sum_track_noi = sum(track_noi,2);
    sum_real_noi = sum(real_noi,2);
    sum_pro = sum_real_noi ./ sum_track_noi;

xout1 = istft(f_sig,FrameLen,0.5 * FrameLen,1);
audiowrite('out1.wav',int16(xout1),fs);
close(hwt)
fclose('all');