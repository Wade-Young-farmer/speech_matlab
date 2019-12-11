function [Y] = wpe_process(X, frame_len, reverb_time)
% X - input subband signal, mic_num * subband_num * frame_num
% frame_len - frame length in milliseconds
% reverb_time - reverberation time (T60) in milliseconds
% Y - out subband signal, mic_num * subband_num * frame_num

mic_num = size(X, 1);
subband_num = size(X, 2);
frame_num = size(X, 3);

% minimum delay parameter (first tap delay)
delta = 20;
min_delay = round(delta/frame_len);

% filter length
% h_size = round((reverb_time / frame_len - 1) / (mic_num - 1));
% h_size_allband = round((reverb_time / frame_len - 1) / (mic_num - 1));
% h_size = [h_size_allband, h_size_allband];
h_size = [7, 7];

% maximum delay parameter
max_delay = h_size + min_delay;

% forgetting factor
gamma = 0.975;
% gamma = 0.98;

% decay rate
T60 = reverb_time;
decay_rate = 3 / (T60 * log10(exp(1)));

% time span outside of which source signals do not have significant
%   auto-correlation coefficients
theta = 50;
theta_frames = round(theta / frame_len);
% beta = exp(1)^(-2*decay_rate*theta_frames);
beta = exp(1)^(-2*decay_rate*theta);

Y = zeros(mic_num, subband_num, frame_num);
Y(:, :, 1:max(max_delay)-1) = X(:, :, 1:max(max_delay)-1);

H = zeros(mic_num * max(h_size), 2, frame_num);
% H = zeros((mic_num-1) * max(h_size), 2, frame_num);
HH = zeros(mic_num * max(h_size), frame_num);

FBF_ONLY_FOR_REF_FLAG = 0;
% RLS style
USE_SFT_RLS = 0;
USE_NORMAL_RLS = 1;

w_debug = zeros(max(h_size), frame_num);

for ss = 1:subband_num
    %% SFT-RLS
    if (USE_SFT_RLS == 1)
        %% SFT-RLS
        H_SIZE = max(h_size);
        X_S = squeeze(X(:, ss, :));
        %初始化数值
        lambda = 0.95;                %定义遗忘因子1-0.4/(N+1)
        delta = 0.01;
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
        
        for ff = max_delay+1:frame_num
            X_SF = X_S(ff);                                                     % 参考信号， (mic_num) * 1
            X_SF_delta = X_S(ff-min_delay:-1:ff-max_delay);
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
            Y(:, ss, ff) = e;
%             if(ss==10)
%                 w_debug(:,ff) = w;
%                 if (ff < 28)
%                     ff
%                     gamma_last_3
%                     xi_b_last
%                     epsilon
%                 end
%             end
        end 
        disp(['process subband ' int2str(ss) '/' int2str(subband_num)]);  
        
    elseif (USE_NORMAL_RLS == 1)
        %% Normal-RLS  
        X_S = squeeze(X(:, ss, :));
        if (ss < ceil(subband_num/2))
            H_SIZE = h_size(1);
            MAX_DELAY = max_delay(1);
        else
            H_SIZE = h_size(2);
            MAX_DELAY = max_delay(2);
        end
        
        if (FBF_ONLY_FOR_REF_FLAG == 1)
            mic_num = 6;
            matrix = 1 * eye(mic_num * H_SIZE, mic_num * H_SIZE);
            h = zeros(mic_num * H_SIZE, mic_num);

            for ff = MAX_DELAY:frame_num
                X_SF = X_S(:, ff);
                X_SF(2) = [];
                X_SF_delta = X_S(:, ff-min_delay:-1:ff-MAX_DELAY+1);
                X_SF_delta(1,:) = [];
                X_SF_delta = reshape(X_SF_delta.', length(X_SF_delta(:)), 1);

                lamda = (1/mic_num) * sum(abs(X_SF).^2);
                k = matrix * X_SF_delta/(gamma * lamda + X_SF_delta' * matrix * X_SF_delta);
                matrix = (1/gamma) * (matrix - k * X_SF_delta' * matrix);

                X_pred = X_SF - h' * X_SF_delta;
                h = h + k * X_pred';
                X_pred = X_SF - h' * X_SF_delta;

                if (sum(abs(X_pred)) > sum(abs(X_SF)))
                    X_pred = X_SF;
                end

                Y(1:6, ss, ff) = X_pred;
%                 Y(:, ss, ff) = X_pred;
                
                % 滤波器收敛分析
                if (ss == 10)
                    H(:, :, ff) = h(:,1:2);
                end
            end       
        else
        %     matrix = eps * eye(mic_num * h_size, mic_num * h_size);
            matrix = 10000 * eye(mic_num * H_SIZE, mic_num * H_SIZE);
            h = zeros(mic_num * H_SIZE, mic_num);

        %     data_buf = zeros(mic_num, 6);
        %     data_buf_id = 1;

            for ff = MAX_DELAY:frame_num
                if (mic_num == 1)
                    X_SF = X_S(ff);
                    X_SF_delta = X_S(ff-min_delay:-1:ff-MAX_DELAY+1);
                    X_SF_delta = reshape(X_SF_delta.', length(X_SF_delta(:)), 1);
                else
                    X_SF = X_S(:, ff);
                    X_SF_delta = X_S(:, ff-min_delay:-1:ff-MAX_DELAY+1);
                    X_SF_delta = reshape(X_SF_delta.', length(X_SF_delta(:)), 1);
                end

        %     	lamda = (1/mic_num) * (sum(abs(X_SF).^2) - beta * sum(abs(X_S(:, ff-theta_frames)).^2));
        %         lamda = abs(lamda);

        %         %初始化权重系数
        %         if(ff <= max_delay + 5)
        %             lamda = (1/(mic_num*6)) * sum(sum(abs(X_S(:, ff:-1:ff-5)).^2));
        %         end

                lamda = (1/mic_num) * sum(abs(X_SF).^2);
                k = matrix * X_SF_delta/(gamma * lamda + X_SF_delta' * matrix * X_SF_delta);
                matrix = (1/gamma) * (matrix - k * X_SF_delta' * matrix);

                X_pred = X_SF - (X_SF_delta' * h).';
                h = h + k * X_pred.';
                X_pred = X_SF - (X_SF_delta' * h).';

                if (sum(abs(X_pred)) > sum(abs(X_SF)))
                    X_pred = X_SF;
                end

                Y(:, ss, ff) = X_pred;

                % 滤波器收敛分析
                if (ss == 10)
                    HH(:, ff) = h(:,1);
                end

        %         %更新权重系数
        %         if (sum(abs(X_pred)) > sum(abs(X_SF)))
        %             data_buf(:, data_buf_id) = X_SF;
        %         else
        %             data_buf(:, data_buf_id) = X_pred;
        %         end
        %         data_buf_id = data_buf_id + 1;
        %         if (data_buf_id == 7) 
        %             data_buf_id = 1;
        %         end
        %         if (ff >=  max_delay + 5)
        %             lamda = (1/(mic_num*6)) * sum(sum(abs(data_buf).^2));
        %         end
            end
        end
        disp(['process subband ' int2str(ss) '/' int2str(subband_num)]);  
    end
end

figure(3);clf
for i=7:7:mic_num*max(h_size)
    plot(abs(HH(i,1:3000)));
    hold on;
end
% save H_low10_high10_6mic.mat HH
end
