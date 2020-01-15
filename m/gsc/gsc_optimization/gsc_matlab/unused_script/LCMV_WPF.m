clear;
clc;

Br=0.99;
Coe=0.2;
% sm_fv=0.9512;
sm_fv=0.9355;

fs=16000;
c = 340;
MicNum = 6;
      
f_len=1024;    %帧长
use_flen=f_len/2+1;
f_step=512;           %帧移
FreMin = floor(0 / fs * f_len) + 1;  %1
FreMax = floor(fs/2 / fs * f_len) + 1;  %513
FreLen = FreMax - FreMin + 1;
hamm_win = hamming(f_len);
fre_disp = (FreMin : FreMax) / f_len * fs;

R_cir = 0.03;
R_tar = 2;
% sig_angle = -100 / 180 * pi;
% dist_angle = -40 / 180 * pi;
sig_angle = -100 / 180 * pi;
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

sig_steering = zeros(MicNum, FreLen);
dist_steering = zeros(MicNum, FreLen);
for f = 1 : FreLen
    fre = (f + FreMin-1) / f_len * fs;
    omiga = 2 * pi * fre;
	sig_steering(:, f) = exp(-j * omiga * sig_mic_tao);
    dist_steering(:, f) = exp(-j * omiga * dist_mic_tao);
end

H1=ones(6,use_flen);
H2=ones(6,use_flen);


fp=fopen('D:\代码\matlab程序\result\superBF.pcm','wb');

% load MicAuData10dBRt03s.mat;          %输入语音为列向量
% x(:,1)=32768*audioread('D:\代码\gsc_matlab\gsc_matlab\xiaoyu_data\rec_mic_0_0.pcm');
% x(:,2)=32768*audioread('D:\代码\gsc_matlab\gsc_matlab\xiaoyu_data\rec_mic_1_0.pcm');
% x(:,3)=32768*audioread('D:\代码\gsc_matlab\gsc_matlab\xiaoyu_data\rec_mic_2_0.pcm');
% x(:,4)=32768*audioread('D:\代码\gsc_matlab\gsc_matlab\xiaoyu_data\rec_mic_3_0.pcm');
% x(:,5)=32768*audioread('D:\代码\gsc_matlab\gsc_matlab\xiaoyu_data\rec_mic_4_0.pcm');
% x(:,6)=32768*audioread('D:\代码\gsc_matlab\gsc_matlab\xiaoyu_data\rec_mic_5_0.pcm');
% filepath = 'D:\代码\gsc_matlab\gsc_matlab\xiaoyu_data\';
filepath = 'D:\代码\gsc_yuzhou\gsc_matlab\xiaoyu_data\';
wavename = 'rec_mic_0_0.pcm';
pathInfo = dir([filepath wavename]);
SigLen = pathInfo.bytes/2;
% SigLen =3034395;
for m = 1 : MicNum
    wavename = ['rec_mic_' num2str(m-1) '_0.pcm'];
    pathInfo = dir([filepath wavename]);
    L = pathInfo.bytes/2;
    if L > SigLen
        L = SigLen;
    end
    fid = fopen([filepath wavename],'rb');
    x(:,m) = fread(fid, [1,L], 'int16');
    fclose(fid);
end
tic
% SigLen = 178 * fs;
frame_num=fix((length(x(:,1))-f_step)/f_step);  

sig_fix = zeros(MicNum, SigLen);
sig_beam = zeros(1, SigLen);
dist_fix = zeros(MicNum, SigLen);
dist_beam = zeros(1, SigLen);

len1024 = 1024;
len513 = 513;

ctrl_sig_beam = zeros(1, len513);
ctrl_sig_fix = zeros(MicNum, len513);
ctrl_dist_beam = zeros(1, len513);
ctrl_dist_fix = zeros(MicNum, len513);

ctrl_snr_disp = zeros(frame_num, len513);
ctrl_nsr_disp = zeros(frame_num, len513);
ctrl_snr_alpha = zeros(1, len513);
ctrl_nsr_alpha = zeros(1, len513);

In_voice1=enframe(x(:,1),sqrt(hanning(f_len,'periodic')),f_step);         %In_voice1~6分别存6路语音，为矩阵，每一行对应一帧语音
In_voice2=enframe(x(:,2),sqrt(hanning(f_len,'periodic')),f_step);
In_voice3=enframe(x(:,3),sqrt(hanning(f_len,'periodic')),f_step);
In_voice4=enframe(x(:,4),sqrt(hanning(f_len,'periodic')),f_step);
In_voice5=enframe(x(:,5),sqrt(hanning(f_len,'periodic')),f_step);
In_voice6=enframe(x(:,6),sqrt(hanning(f_len,'periodic')),f_step);
    
    
   
Last_voice=0;    %保存上一帧后512个数据
sig_index=ones(1,use_flen);
noi_index=ones(1,use_flen);
for i=1:20
    fv(:,:,i)=eye(6,use_flen);
end
sig_buffer=zeros(6,use_flen,20);
% faiv=zeros(6,6,use_flen);
for i=1:use_flen
    faiv_inv(:,:,i)=eye(6,6);
end
hwt = waitbar(0, 'GSC process'); 
for frame=3000:frame_num-5000
    waitbar(frame/frame_num, hwt);
    
    % data prepare
    f_voice1=fft(In_voice1(frame,:));        %f_voice1~6中存每一路第frame帧傅里叶变换的结果
    f_voice2=fft(In_voice2(frame,:));        %为行向量
    f_voice3=fft(In_voice3(frame,:)); 
    f_voice4=fft(In_voice4(frame,:)); 
    f_voice5=fft(In_voice5(frame,:)); 
    f_voice6=fft(In_voice6(frame,:));
    ll = (frame-1)*f_step + 1 : (frame-1)*f_step + f_len;
    frame_data = x(ll,:);
%     frame_data = frame_data';
    for m = 1 : MicNum
        frame_data(:,m) = frame_data(:,m) .* hamm_win;
    end
    fft_data = fft(frame_data);
    
        % signal steering and fix beam
    for m = 1 : MicNum
        fft_fix = zeros(f_len, 1);
        fft_fix(1 : FreLen) = fft_data(1 : FreLen, m) .* sig_steering(m, 1 : FreLen)';
        fft_fix(f_len/2+2:f_len) = conj(fft_fix(f_len/2:-1:2));
        fix_temp = real(ifft(fft_fix));
        fix_temp = fix_temp';
        sig_fix(m, ll) =fix_temp;
    end
    sig_beam(ll) = mean(sig_fix(:, ll));
    
        % disturb steering and fix beam
    for m = 1 : MicNum
        fft_fix = zeros(f_len, 1);
        fft_fix(1 : FreLen) = fft_data(1 : FreLen, m) .* dist_steering(m, 1 : FreLen)';
        fft_fix(f_len/2+2:f_len) = conj(fft_fix(f_len/2:-1:2));
        fix_temp = real(ifft(fft_fix));
        fix_temp = fix_temp';
        dist_fix(m, ll) =  fix_temp;
    end
    dist_beam(ll) = mean(dist_fix(:, ll));
    
        ctrl_sig_beam_fft = fft(sig_beam(ll));
        ctrl_sig_fix_fft = fft(sig_fix(:, ll), [], 2);
        ctrl_dist_beam_fft = fft(dist_beam(ll));
        ctrl_dist_fix_fft = fft(dist_fix(:, ll), [], 2);
        
        ctrl_sig_null_eng = zeros(1, len513);
        ctrl_dist_null_eng = zeros(1, len513);
        for m = 1 : MicNum
            temp = ctrl_sig_beam_fft(1 : len513) - ctrl_sig_fix_fft(m, 1 : len513);
            ctrl_sig_null_eng = ctrl_sig_null_eng + abs(temp).^2 / MicNum;
            temp = ctrl_dist_beam_fft(1 : len513) - ctrl_dist_fix_fft(m, 1 : len513);
            ctrl_dist_null_eng = ctrl_dist_null_eng + abs(temp).^2 / MicNum;
        end
        for n = 1 : len513
            if ctrl_sig_null_eng(n) < 0.001
                ctrl_sig_null_eng(n) = 0.001;
            end
            if ctrl_dist_null_eng(n) < 0.001
                ctrl_dist_null_eng(n) = 0.001;
            end
        end
        ctrl_snr = ctrl_dist_null_eng ./ ctrl_sig_null_eng;
        ctrl_nsr = ctrl_sig_null_eng ./ ctrl_dist_null_eng;
        for n = 8 : len513
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
            snr_temp(n) = ctrl_snr_alpha(n) / (ctrl_snr_alpha(n) + ctrl_nsr_alpha(n));
            
        end
        for n = 1 : 7
            ctrl_snr_alpha(n)=ctrl_snr_alpha(n+7);
            ctrl_nsr_alpha(n)=ctrl_nsr_alpha(n+7);
        end
        ctrl_snr_disp(frame, :) = ctrl_snr_alpha;
        ctrl_nsr_disp(frame, :) = ctrl_nsr_alpha;
        
        test_snr=ctrl_snr_alpha(1:use_flen)./ctrl_nsr_alpha(1:use_flen);
        test_sum_snr=sum(test_snr);
        
        
        for i=1:use_flen
            if(test_snr(i)<0.8)  %as noise
                tmp_fv=[f_voice1(i);f_voice2(i);f_voice3(i);f_voice4(i);f_voice5(i);f_voice6(i)];
                if(noi_index(i)<=20)
                    fv(:,i,noi_index(i))=tmp_fv;
                    noi_index(i)=noi_index(i)+1;
                else
                    for j=1:19
                        fv(:,i,j)=fv(:,i,j+1);
                    end
                    fv(:,i,20)=tmp_fv;
                    noi_index(i)=21;
                end
                
                tmp_faivinv=faiv_inv(:,:,i)/sm_fv;
                faiv_inv(:,:,i)=tmp_faivinv-(1-sm_fv)*(tmp_faivinv*(tmp_fv*tmp_fv')*tmp_faivinv)/(1+(1-sm_fv)*tmp_fv'*tmp_faivinv*tmp_fv);
            end 
            
            if(test_snr(i)>1.2)  %as signal
                tmp_sig=[f_voice1(i);f_voice2(i);f_voice3(i);f_voice4(i);f_voice5(i);f_voice6(i)];
                if(sig_index(i)<=20)
                    sig_buffer(:,i,sig_index(i))=tmp_sig;
                    sig_index(i)=sig_index(i)+1;
                else
                    for j=1:19
                        sig_buffer(:,i,j)=sig_buffer(:,i,j+1);
                    end
                    sig_buffer(:,i,20)=tmp_sig;
                    sig_index(i)=21;
                end
            end
        end
            
                       
    %***************************计算H1***************************************
        fai1=zeros(6,use_flen);      %变量初始化
        fai2=zeros(1,use_flen); 
        fai3=zeros(6,use_flen); 
        fai4=zeros(1,use_flen);   
        fai5=zeros(1,use_flen); 
        
    for frm=1:use_flen
         if(sig_index(frm)>2)
            for i=1:sig_index(frm)-1   %控制帧循环

            fai2(frm)=fai2(frm)+sig_buffer(1,frm,i)*conj(sig_buffer(1,frm,i));
            fai4(frm)=fai4(frm)+sig_buffer(1,frm,i)*conj(sig_buffer(1,frm,i))*sig_buffer(1,frm,i)*conj(sig_buffer(1,frm,i));

                for m=1:6
                    fai1(m,frm)=fai1(m,frm)+sig_buffer(1,frm,i)*conj(sig_buffer(1,frm,i))*sig_buffer(m,frm,i)*conj(sig_buffer(1,frm,i));
                    fai3(m,frm)=fai3(m,frm)+sig_buffer(m,frm,i)*conj(sig_buffer(1,frm,i));
                end
            end
            fai1(:,frm)=fai1(:,frm)/(sig_index((frm))-1);
            fai2(frm)=fai2(frm)/(sig_index((frm))-1);
            fai3(:,frm)=fai3(:,frm)/(sig_index((frm))-1);
            fai4(frm)=fai4(frm)/(sig_index((frm))-1);
            fai5(frm)=fai2(frm)*conj(fai2(frm));

            for m=1:6
                H1(m,frm)=(fai1(m,frm)-fai2(frm)*fai3(m,frm))/(fai4(frm)-fai5(frm));
            end

         end
    end
    %******************************************************************
         
    
    %***************************计算H2***************************************
        fai1=zeros(6,use_flen);      %变量初始化
        fai2=zeros(1,use_flen); 
        fai3=zeros(6,use_flen); 
        fai4=zeros(1,use_flen); 
        fai5=zeros(1,use_flen); 
        
    for frm=1:use_flen
        if(noi_index(frm)>2)
            for i=1:noi_index(frm)-1   %控制帧循环

            fai2(frm)=fai2(frm)+fv(1,frm,i)*conj(fv(1,frm,i));
            fai4(frm)=fai4(frm)+fv(1,frm,i)*conj(fv(1,frm,i))*fv(1,frm,i)*conj(fv(1,frm,i));

                for m=1:6
                    fai1(m,frm)=fai1(m,frm)+fv(1,frm,i)*conj(fv(1,frm,i))*fv(m,frm,i)*conj(fv(1,frm,i));
                    fai3(m,frm)=fai3(m,frm)+fv(m,frm,i)*conj(fv(1,frm,i));
                end
            end
            fai1(:,frm)=fai1(:,frm)/(noi_index(frm)-1);
            fai2(frm)=fai2(frm)/(noi_index(frm)-1);
            fai3(:,frm)=fai3(:,frm)/(noi_index(frm)-1);
            fai4(frm)=fai4(frm)/(noi_index(frm)-1);
            fai5(frm)=fai2(frm)*conj(fai2(frm));
            for m=1:6
                H2(m,frm)=(fai1(m,frm)-fai2(frm)*fai3(m,frm))/(fai4(frm)-fai5(frm));
            end
        end
    end
    %******************************************************************

    
        for frm=1:use_flen       %代表频率值
        %     B(:,:,fre)=zeros(6,4);
            if((noi_index(frm)>2)&&(sig_index(frm)>2))
                G(:,:,frm)=[H1(:,frm) H2(:,frm)];
            %     H0(:,:,fre)=G(:,:,fre)*inv(G(:,:,fre)'*G(:,:,fre));
            %     B(1:2,1,fre)=-inv(G(1:2,1:2,fre))'*G(3,:,fre)';
            %     B(1:2,2,fre)=-inv(G(1:2,1:2,fre))'*G(4,:,fre)';
            %     B(1:2,3,fre)=-inv(G(1:2,1:2,fre))'*G(5,:,fre)';
            %     B(1:2,4,fre)=-inv(G(1:2,1:2,fre))'*G(6,:,fre)';
            %     B(3:6,1:4,fre)=ones(4,4);
            %     B_F=norm(B(:,:,fre),2);
            %     B(:,:,fre)=B(:,:,fre)/B_F;
            %     matBB=B(:,:,fre)'*faiv(:,:,fre)*B(:,:,fre);
            %     eB=max(eig(matBB))/10000;
            %     BH(:,:,fre)=B(:,:,fre)*inv(matBB+eB*ones(4,4))*B(:,:,fre)'*faiv(:,:,fre)*H0(:,:,fre);
            %     Hlcmv(:,:,fre)=H0(:,:,fre)-BH(:,:,fre);
                faivre(:,:,frm)=inv(G(:,:,frm)'*faiv_inv(:,:,frm)*G(:,:,frm));
                Hlcmv(:,:,frm)=faiv_inv(:,:,frm)*G(:,:,frm)*faivre(:,:,frm);

                firH(:,frm)=Hlcmv(:,1,frm);
                secH(:,frm)=Hlcmv(:,2,frm);
            else
                firH(:,frm)=ones(6,1);
                secH(:,frm)=ones(6,1);
                faivre(:,:,frm)=eye(2,2);
            end
        end

    
    fai_sf=zeros(2,2,use_flen);
    
    
    F_voice1=f_voice1(1:use_flen).*conj(firH(1,1:use_flen));     %第一个说话人
    F_voice2=f_voice2(1:use_flen).*conj(firH(2,1:use_flen));
    F_voice3=f_voice3(1:use_flen).*conj(firH(3,1:use_flen));
    F_voice4=f_voice4(1:use_flen).*conj(firH(4,1:use_flen));
    F_voice5=f_voice5(1:use_flen).*conj(firH(5,1:use_flen));
    F_voice6=f_voice6(1:use_flen).*conj(firH(6,1:use_flen));
    
    S_voice1=f_voice1(1:use_flen).*conj(secH(1,1:use_flen));      %第二个说话人
    S_voice2=f_voice2(1:use_flen).*conj(secH(2,1:use_flen));
    S_voice3=f_voice3(1:use_flen).*conj(secH(3,1:use_flen));
    S_voice4=f_voice4(1:use_flen).*conj(secH(4,1:use_flen));
    S_voice5=f_voice5(1:use_flen).*conj(secH(5,1:use_flen));
    S_voice6=f_voice6(1:use_flen).*conj(secH(6,1:use_flen));
      
    Sy_Fvoice=(F_voice1+F_voice2+F_voice3+F_voice4+F_voice5+F_voice6)/6;       %暂算513点,行向量
    Sy_Svoice=(S_voice1+S_voice2+S_voice3+S_voice4+S_voice5+S_voice6)/6;
    
    Slcmv=[Sy_Fvoice;Sy_Svoice];         %两行分别存两路语音
    
    
    if(frame==3000)
        Wp_voice=Slcmv;  %zeros(2,513);   
    end
  

    for fre=1:use_flen
        fai_sf(1,1,fre)=Br*(Wp_voice(1,fre)*conj(Wp_voice(1,fre)))+(1-Br)*max(Sy_Fvoice(fre)*conj(Sy_Fvoice(fre))-real(faivre(1,1,fre)),0);
        fai_sf(2,2,fre)=Br*(Wp_voice(2,fre)*conj(Wp_voice(2,fre)))+(1-Br)*max(Sy_Svoice(fre)*conj(Sy_Svoice(fre))-real(faivre(2,2,fre)),0);
        sfai=fai_sf(:,:,fre)+faivre(:,:,fre);

%         r=max(abs(eig(sfai)));
        Hwpf(:,:,fre)=(1-Coe)*fai_sf(:,:,fre)*inv(sfai)+Coe*eye(2,2);
        Wp_voice(:,fre)=Hwpf(:,:,fre)*Slcmv(:,fre);

    end

    W_voice=Wp_voice(1,:);
    W_voice(1)=real(W_voice(1));
    W_voice(use_flen)=real(W_voice(use_flen));
    W_voice(f_len/2+2:f_len)=conj(fliplr(W_voice(2:f_len/2)));
    
    T_voice=real(ifft(W_voice)).*sqrt(hanning(f_len,'periodic')).';
   
    Overlap_fvoice=0;
    Overlap_voice=Last_voice+T_voice(1:f_len/2);   
    fwrite(fp, 20*Overlap_voice,'int16');
%     fprintf(fp,'%d\r\n',20*Overlap_voice);
%     audiowrite('D:\代码\matlab程序\superBF.wav', Overlap_voice, fs);
        
    Last_voice=T_voice(use_flen:f_len);           %保存上一帧后512个数据    
end
close(hwt)
toc
fclose(fp);
fclose('all');