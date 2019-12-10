%FTF快速横向滤波算法
clc;
clear all;
close all;
%************************生成仿真信号**************************************
Fs = 10000;                                                     %设置采样频率
t = 0:1/Fs:3.5;  
t = t';
Size_t = size(t,1);
F1 = 2;
F2 = 10;
F3 = 20;
F4 = 1000;
Signal = sin(2*pi*F1*t) + 0.5*sin(2*pi*F2*t) + 0.25*sin(2*pi*F3*t); %生成信号
noise_amp = 1;                                           %定义噪声的标准差
noise1 = noise_amp*randn(Size_t,1);                        %生成高斯白噪声
noise2 = noise_amp*randn(Size_t,1);
noise3 = 5*sin(2*pi*F4*t);

noise = noise2;
Signal_noise = Signal + 0.2*noise;                           %加入高斯白噪声
Signal_noise(2:end) = Signal_noise(2:end) + 0.15*noise(1:end-1);
Signal_noise(3:end) = Signal_noise(3:end) + 0.1*noise(1:end-2);

subplot(2,1,1);
plot(t,Signal);
title('原始信号');
subplot(2,1,2);
plot(t,Signal_noise);
title('加入干扰噪声的信号');

%% 子带分析
addpath('./filter_bank');
frame_len = 10;
frame_size = Fs * frame_len / 1000;
over_sample_ratio = 2;
subband_num = frame_size * over_sample_ratio / 2 + 1;
% frame_num = floor(size(x, 1) / frame_size);
frame_num = ceil(size(Signal_noise, 1) / frame_size);
delay = 1;
h_size = frame_size * over_sample_ratio * 3;
h = filter_bank_win(h_size, delay);

% 参考信号
mic_num = 1;
X = zeros(mic_num, subband_num, frame_num);
for ii = 1:mic_num
    X(ii, :, :) = filter_bank_analyze(Signal_noise(:, ii), frame_size, over_sample_ratio, h);
end
XX = X;
% 信号
Z = zeros(mic_num, subband_num, frame_num);
for ii = 1:mic_num
    Z(ii, :, :) = filter_bank_analyze(noise(:, ii), frame_size, over_sample_ratio, h);
end
ZZ = Z;

%% RLS处理
% tic
% Y = stf_rls_complex_monoral(XX, frame_len);
% toc

H_SIZE = 3;
min_delay = 0;
max_delay = H_SIZE;
Eta_out = zeros(subband_num, frame_len);
for ss=1:subband_num
        X_S = squeeze(X(:, ss, :));
        Z_S = squeeze(Z(:, ss, :));
        %初始化数值
        lambda = 0.9;%1.0001;                %定义遗忘因子1-0.4/(N+1)
        delta = 0.0000001;%0.01;
        w_f_last = zeros(H_SIZE*mic_num,mic_num);
        w_b_last = zeros(H_SIZE*mic_num,mic_num);
        w_last = zeros(H_SIZE*mic_num,mic_num);
        Phi_last = zeros(H_SIZE*mic_num,mic_num);
        gamma_last_3 = ones(mic_num,1);
        xi_f_last = delta * ones(mic_num,1);
        xi_b_last = delta * ones(mic_num,1);
        k1 = 1.5;
        k2 = 2.5;
        k3 = 1;
        
        for ff = H_SIZE+1:frame_num
            X_SF = X_S(ff);                                                     % 参考信号， (mic_num) * 1
            X_SF_delta = Z_S(ff-min_delay:-1:ff-max_delay);
            X_SF_delta = reshape(X_SF_delta.', length(X_SF_delta(:)), 1);          % 历史信号， (mic_num*order) * 1
         
            %算法主体
            x_N_1 = X_SF_delta;
            d = X_SF;
            
            e_f = [1, -1*w_f_last'] * x_N_1;    % e_f = x_N_1' * [1;-w_f_last];                                            %(1)(1)
            epsilon_f = e_f * gamma_last_3;                                          %(2)(2)
            Phi_N_1 = [zeros(mic_num,mic_num);Phi_last] + e_f/(lambda * xi_f_last) * [ones(mic_num,mic_num);-w_f_last];    % Phi_N_1 = [0;Phi_last] + e_f/(lambda * xi_f_last)*[1;-w_f_last];         %(3)(5)
            Phi_N_1_0 = Phi_N_1(1:mic_num);
            Phi_N_1_N_1 = Phi_N_1(end-mic_num+1:end);
            gamma_N_1_1 = 1/(1/gamma_last_3 + Phi_N_1_0 * e_f);                      %(4)
            xi_f = 1/(1/(lambda * xi_f_last) - (Phi_N_1_0 * Phi_N_1_0) * gamma_N_1_1); %(5)
            w_f = w_f_last + Phi_last * epsilon_f;                                 %(6)(4)
%             w_f = w_f_last + Phi_last .* repmat(epsilon_f', h_size*mic_num, 1);                                    %(6)(4)
            e_b_1 = lambda * Phi_N_1_N_1 * xi_b_last;                                %(7)(7)
            e_b_2 = [-w_b_last',ones(mic_num,mic_num)] * x_N_1;                                          %(8)
            e_b_3_1 = e_b_2 * k1 + e_b_1 * (1-k1);                                   %(9)
            e_b_3_2 = e_b_2 * k2 + e_b_1 * (1-k2);
            e_b_3_3 = e_b_2 * k3 + e_b_1 * (1-k3);
            gamma_2 = 1/(1/gamma_N_1_1 - Phi_N_1_N_1 * e_b_3_3);                     %(10)
            epsilon_b_3_1 = e_b_3_1 * gamma_2;                                       %(11)
            epsilon_b_3_2 = e_b_3_2 * gamma_2;
            xi_b = lambda * xi_b_last + epsilon_b_3_2 * e_b_3_2;                     %(12)(10)
            Phi = Phi_N_1 - [-w_b_last;ones(mic_num,mic_num)] * Phi_N_1_N_1;                             %(13)(11)
            Phi = Phi(1:end-mic_num); 
            w_b = w_b_last + Phi * epsilon_b_3_1;                                    %(14)(12)
            x = x_N_1(1:H_SIZE*mic_num);
            gamma_3 = 1/(ones(mic_num,1) + Phi' * x);                                              %(15)
            %联合过程估计
            y = w_last'* x;                                                  
            e = d - y;                                                               %(16)(13)
            epsilon = e * gamma_3;                                                   %(17)(14)
            w = w_last + Phi * epsilon;                                              %(18)(15)
            %变量更替
            xi_f_last = xi_f;
            w_f_last = w_f;
            gamma_last_3 = gamma_3;
            xi_b_last = xi_b;
            Phi_last = Phi;
            w_b_last = w_b;
            w_last = w;
            %滤波结果存储
            Eta_out(ss, ff) = e;
        end 
        disp(['process subband ' int2str(ss) '/' int2str(subband_num)]);  
end

for jj = 1:mic_num
    YY = Eta_out;
    y = filter_bank_synthesis(YY, frame_size, over_sample_ratio, h);
    figure(jj+1);
    plot(y);
    title('输出误差');
%     out_filename = ['rls_out' '_0301_' int2str(jj+2) '.wav'];
%     out_file = fullfile(filepath, out_filename);
%     y = y ./ 32768;
%     audiowrite(out_file, y, fs);
end

   



