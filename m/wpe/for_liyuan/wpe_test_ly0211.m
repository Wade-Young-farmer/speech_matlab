clc;
clear all;
close all;

if (1)
    mic_num = 1;

    filepath = 'E:\Baidu_Work\项目\听清\识别在线测试\20181220\doris_rawcut_1220\xf_doris对齐\doris_rawcut\medium\LQ_LZ_ask3_M\';

    if (mic_num == 1)
%         waveName = 'recog_0218_afterFBF_noAGC_forWPE.pcm';
        waveName = 'recog_0222_LOC_FBF.pcm';
        pathInfo = dir([filepath waveName]);
        L = pathInfo.bytes/2;
        fid = fopen([filepath waveName],'rb');
        sigMic(:,1) = fread(fid,[L,1],'int16');
        fclose(fid);
    elseif (mic_num == 2)
        waveName = 'orig3.pcm';
        pathInfo = dir([filepath waveName]);
        L = pathInfo.bytes/2;
        fid = fopen([filepath waveName],'rb');
        sigMic(:,1) = fread(fid,[L,1],'int16');
        fclose(fid);
        waveName = 'orig6.pcm';
        pathInfo = dir([filepath waveName]);
        L = pathInfo.bytes/2;
        fid = fopen([filepath waveName],'rb');
        sigMic(:,2) = fread(fid,[L,1],'int16');
        fclose(fid);
    elseif (mic_num == 3)
        waveName = 'orig3.pcm';
    %     waveName = 'orig4.pcm';
        pathInfo = dir([filepath waveName]);
        L = pathInfo.bytes/2;
        fid = fopen([filepath waveName],'rb');
        sigMic(:,1) = fread(fid,[L,1],'int16');
        fclose(fid);
        waveName = 'orig5.pcm';
    % waveName = 'orig6.pcm';
        pathInfo = dir([filepath waveName]);
        L = pathInfo.bytes/2;
        fid = fopen([filepath waveName],'rb');
        sigMic(:,2) = fread(fid,[L,1],'int16');
        fclose(fid);
        waveName = 'orig7.pcm';
    %     waveName = 'orig8.pcm';
        pathInfo = dir([filepath waveName]);
        L = pathInfo.bytes/2;
        fid = fopen([filepath waveName],'rb');
        sigMic(:,3) = fread(fid,[L,1],'int16');
        fclose(fid);
    elseif (mic_num == 6)
        for nn = 1:mic_num
            waveName = ['orig' int2str(nn+2) '.pcm'];
            pathInfo = dir([filepath waveName]);
            L = pathInfo.bytes/2;
            fid = fopen([filepath waveName],'rb');
            sigMic(:,nn) = fread(fid,[L,1],'int16');
            fclose(fid);
        end
    elseif (mic_num == 7)
%         waveName = 'recog_0220_onlyFBF.pcm';
        waveName = 'recog_0222_LOC_FBF.pcm';
        pathInfo = dir([filepath waveName]);
        L1 = pathInfo.bytes/2;
        fid = fopen([filepath waveName],'rb');
        sigMic(:,1) = fread(fid,[L1,1],'int16');
        fclose(fid);
        for nn = 1:mic_num-1
            waveName = ['orig' int2str(nn+2) '.pcm'];
            pathInfo = dir([filepath waveName]);
            L = pathInfo.bytes/2;
            L_min = min(L1,L);
            fid = fopen([filepath waveName],'rb');
            sigMic(L_min+1:end, :) = [];
            sigMic(:,nn+1) = fread(fid,[L_min,1],'int16');
            fclose(fid);
        end
    end

    fs = 16000;

    dbstop if error;

    addpath('./filter_bank');

    % file_path = './wav_sample';
    % in_file = fullfile(file_path, 'sample_4ch.wav');
    reverb_time = 500;

    % [x, fs] = audioread(in_file);
    % mic_num = size(x, 2);
    x = sigMic;

    frame_len = 10;
    frame_size = fs * frame_len / 1000;
    over_sample_ratio = 2;
    subband_num = frame_size * over_sample_ratio / 2 + 1;
    % frame_num = floor(size(x, 1) / frame_size);
    frame_num = ceil(size(x, 1) / frame_size);

    delay = 1;
    h_size = frame_size * over_sample_ratio * 3;
    h = filter_bank_win(h_size, delay);

    X = zeros(mic_num, subband_num, frame_num);

    for ii = 1:mic_num
        X(ii, :, :) = filter_bank_analyze(x(:, ii), frame_size, over_sample_ratio, h);
    end

    % % 跑部分数据
    XX = X(:,:,1:6000);
%     XX = X;

%     save XX.mat XX

end

% load XX.mat
% frame_len = 10;
% reverb_time = 500;

tic
Y = wpe_process(XX, frame_len, reverb_time);
toc

% save Y.mat Y

for jj = 1:mic_num
    YY = squeeze(Y(jj, :, :));
    y = filter_bank_synthesis(YY, frame_size, over_sample_ratio, h);
    out_filename = ['wpe_out_rt' int2str(reverb_time) '_0225_' int2str(jj+2) '.wav'];
    out_file = fullfile(filepath, out_filename);
    y = y ./ 32768;
    audiowrite(out_file, y, fs);
end
    
% figure(1);clf
% plot((0:length(x)-1)/fs, x(:, 1), 'b');
% hold on
% plot((0:length(y)-1)/fs, y, 'r');
% grid on



%% 收敛性分析
% load H_low10_high10_6mic.mat
% h_channel1 = imag(squeeze(H(1,1,:)));
% h_channel2 = imag(squeeze(H(1,2,:)));

