function total_input_f_pcm=aec(mic_name, mic_name_2, ref_name, filter_coeff, WEB_RTC_AEC_NL_WEIGHT_CURVE, DOUBLETALK_BAND_TABLE, BAND_TABLE)
[r_enframe, r_f]=hpf(ref_name, filter_coeff);
[x_enframe, x_f]=hpf(mic_name, filter_coeff);
[x2_enframe, x2_f]=hpf(mic_name_2, filter_coeff);

ref_peak=0;
peak_hold_frame=0;

detect_far_end_buffer=zeros(1, 6);
far_end_hold_time=1;
no_ref_count=0;

energy_group=zeros(1,4);
energy_group_peak=zeros(1,4);

stack_ref_low=zeros(32, 20);
stack_ref_hi=zeros(96, 16);

AEC_HPF_COEFF=[0.0, 0.0396, 0.7000, 0.9204];
filter_freeze=0;
freeze_count=0;

fir_coeff_low=zeros(32,20);
adf_coeff_low=zeros(32,20);
fir_coeff_hi=zeros(96, 16);
adf_coeff_hi=zeros(96, 16);
mse_fir=zeros(1, 129);
mse_adf=zeros(1, 129);
mse_mic=zeros(1, 129);

erl_ratio=ones(1,4)*0.25;
dt_flag=0;
fir_update_flag=zeros(1, 129);
fir_out=zeros(1, 129);
fir_out2=zeros(1, 129);
power_ori_mic=zeros(1, 129);
fir_est_ref=zeros(1, 129);
power_echo=zeros(1,129);
% used in fir part, strang meaning
ref_cnt=0;
erl_peak=zeros(1,4);

nlp_parameters = nlp_initial();
nlp_parameters2 = nlp_initial();

dt_spk_peak=zeros(1, 37);
dt_frame_count = 0;
mic_ratio_smooth = 0;
dt_count=0;

total_input_f=zeros(1, 258);
stack_est_ref=zeros(129,6);
dt_hold_age=40;
post_over_drive_sm=0;

fir_coeff_2=zeros(129,6);
adf_coeff_2=zeros(129,6);
mse_fir2=zeros(1, 129);
mse_adf2=zeros(1, 129);
mse_mic2=zeros(1, 129);

spk_t_ns_parameters = noise_initial(16, 100, 200);
for i=1:4
    energy_group_parameters(i)=noise_initial(0.01, 16, 200);
end

for i=1:129
    noise_est_mic_parameters(i)=noise_initial(0.1, 100, 200);
end
doubletalk_1_parameters=noise_initial(0.1, 100, 62);
doubletalk_2_parameters=noise_initial(0.01, 4, 200);
post_parameters=noise_initial(0.01, 1, 100);

hh=waitbar(0, 'data is being processed');
pcm_size=size(x_enframe,1);

total_input_f_pcm=zeros(pcm_size, 258);
for i=1:size(x_enframe,1)
    str=[num2str(i), ' / ', num2str(pcm_size), ' processed']; 
    waitbar(i/pcm_size, hh, str); 
    
    input_t=x_enframe(i,:);
    input_f=x_f(i,:);
    input_f_bak=input_f;

    ref_f=r_f(i,:);
    input_r=ref_f;
    
    % Only apply to single-channel reference signal
    ref_t=r_enframe(i,:);
    if max(abs(ref_t)) > ref_peak
        ref_peak = max(abs(ref_t));
        peak_hold_frame=0;
    else
        peak_hold_frame=peak_hold_frame+1;
        if peak_hold_frame > 100
            ref_peak=(1-0.0952)*ref_peak + 0.0952*max(abs(ref_t));
        end
    end
    
    % Far end talk detect
    magnitude_tmp=mean(abs(ref_t));
    if i == 1
        detect_far_end_buffer = detect_far_end_buffer + magnitude_tmp;
    else
        detect_far_end_buffer(1:5)=detect_far_end_buffer(2:6);
        detect_far_end_buffer(6)=magnitude_tmp;
    end
    
    [spk_t_ns_parameters] = aec_noise_estimation(detect_far_end_buffer(1), spk_t_ns_parameters);
    A=spk_t_ns_parameters.noise_level(1);
    if detect_far_end_buffer(1) > max(10, 2*A)
        far_end_hold_time = 20;
    else
        far_end_hold_time = max(0, far_end_hold_time - 1);
    end

    if far_end_hold_time == 0 && ref_peak < 20
        no_ref_count=min(1000, no_ref_count+1);
    else
        no_ref_count=0;
    end

    ref_e = abs(ref_f).^2;
    for j=1:4
        energy_group(j) = sum(ref_e(BAND_TABLE(2*j-1):BAND_TABLE(2*j)));
        if energy_group(j) > energy_group_peak(j)
            energy_group_peak(j) = energy_group(j);
        else
            energy_group_peak(j) = 0.9048 * energy_group_peak(j) + (1-0.9048) * energy_group(j);
        end
        [energy_group_parameters(j)] = aec_noise_estimation(energy_group(j), energy_group_parameters(j));
    end
    
    stack_ref_low(:, 2:20)=stack_ref_low(:, 1:19);
    stack_ref_low(:, 1)=ref_f(1:32).';
    stack_ref_hi(:,2:16)=stack_ref_hi(:,1:15);
    stack_ref_hi(:,1)=ref_f(33:128).';
    if no_ref_count < 500
        input_f(1:4)=input_f(1:4).*AEC_HPF_COEFF;
        mic_peak=max(abs(input_t));
        % naming not reasonable
        if far_end_hold_time < 20
            ref_cnt=20; 
        else
            ref_cnt=max(0, ref_cnt-1); % there is a far end lasting for 20 frames, would make a zero
        end
        
        if mic_peak > 28000
            filter_freeze=1;
            freeze_count=30;
        else
            freeze_count=max(0, freeze_count-1);
            if freeze_count == 0
                filter_freeze=0;
            end
        end

        for k=1:32
            est_fir_early=stack_ref_low(k,1:4) * fir_coeff_low(k,1:4)';
            est_adf_early=stack_ref_low(k,1:4) * adf_coeff_low(k,1:4)';
            est_fir=stack_ref_low(k,:) * fir_coeff_low(k,:)'; % Matrix Mul, warning of not a scalar
            est_adf=stack_ref_low(k,:) * adf_coeff_low(k,:)';
            est_fir_e = abs(est_fir).^2;
            est_adf_e = abs(est_adf).^2;
            err_fir=input_f(k)-est_fir;
            err_adf=input_f(k)-est_adf;
            input_power=sum(abs(stack_ref_low(k,:)).^2); % sum of buffer of length of 20
            norm_step=0.5/(input_power+0.01);
            err_fir_e=abs(err_fir)^2;
            err_adf_e=abs(err_adf)^2;
            mse_fir(k)=0.9735*mse_fir(k)+(1-0.9735)*err_fir_e;
            mse_adf(k)=0.9735*mse_adf(k)+(1-0.9735)*err_adf_e;
            mse_mic(k)=0.9735*mse_mic(k)+(1-0.9735)*(abs(input_f(k)).^2);
            
            if mse_adf(k) > mse_mic(k) * 8
                adf_coeff_low(k,:)=0;
                mse_adf(k)=mse_mic(k);
                err_adf=input_f(k);
            else
                if mse_adf(k) < min(0.125 * mse_mic(k), 0.5 * mse_fir(k)) && ref_cnt == 0
                    fir_coeff_low(k,:)=adf_coeff_low(k,:);
                    mse_fir(k)=mse_adf(k);
                    err_fir=err_adf;
                end
            end
            
            if mse_fir(k) > mse_mic(k) * 8 && ref_cnt == 0
                fir_coeff_low(k,:)=0;
                mse_fir(k)=mse_mic(k);
                err_fir=input_f(k);
            else
                if mse_fir(k) < min(0.125 * mse_mic(k), 0.5 * mse_adf(k))
                    adf_coeff_low(k,:)=fir_coeff_low(k,:);
                    mse_adf(k)=mse_fir(k);
                    err_adf=err_fir;
                end
            end
            
            if k <= BAND_TABLE(2)
                subband_index=1;
            elseif k <= BAND_TABLE(4)
                subband_index=2;
            elseif k <= BAND_TABLE(6)
                subband_index=3;
            end
            
            if dt_flag == 2
                tmp = 80;
            else
                tmp = 20;
            end
            
            B=noise_est_mic_parameters(k).noise_level(1);
            if input_power * erl_ratio(subband_index) > tmp * B && filter_freeze == 0 % tmp is for dt state threshold
                adf_coeff_low(k,:) = adf_coeff_low(k,:) + norm_step * err_adf' * stack_ref_low(k,:);
                fir_update_flag(k)=1;
            else
                fir_update_flag(k)=0;
            end
            
            adf_early_tmp=input_f(k)-est_adf_early;
            fir_early_tmp=input_f(k)-est_fir_early;
            
            if abs(fir_early_tmp)^2 > abs(adf_early_tmp)^2
                min_power_early=abs(adf_early_tmp)^2;
                early_out=adf_early_tmp;
            else
                min_power_early=abs(fir_early_tmp)^2;
                early_out=fir_early_tmp;
            end

            if err_fir_e > err_adf_e
                power_ori_mic(k)=err_adf_e;
                power_echo(k)=est_adf_e;
                if err_adf_e > min_power_early
                    fir_out(k)=early_out;
                else
                    fir_out(k)=err_adf;
                end 
            else
                power_ori_mic(k)=err_fir_e;
                power_echo(k)=est_fir_e;
                if err_fir_e > min_power_early
                    fir_out(k)=early_out;
                else
                    fir_out(k)=err_fir;
                end
            end
            fir_est_ref(k)=input_f(k)-fir_out(k);
        end
        
        for k = 1:2
            fir_out(k)=0;
        end
          
        % High channel copied from low
        for k=33:128
            est_fir_early=stack_ref_hi(k-32,1:4) * fir_coeff_hi(k-32,1:4)';
            est_adf_early=stack_ref_hi(k-32,1:4) * adf_coeff_hi(k-32,1:4)';
            est_fir=stack_ref_hi(k-32,:) * fir_coeff_hi(k-32,:)'; % Matrix Mul, warning of not a scalar
            est_adf=stack_ref_hi(k-32,:) * adf_coeff_hi(k-32,:)';
            est_fir_e = abs(est_fir).^2;
            est_adf_e = abs(est_adf).^2;
            err_fir=input_f(k)-est_fir;
            err_adf=input_f(k)-est_adf;
            input_power=sum(abs(stack_ref_hi(k-32,:)).^2);
            norm_step=0.5/(input_power+0.01);
            err_fir_e=abs(err_fir)^2;
            err_adf_e=abs(err_adf)^2;
            mse_fir(k)=0.9735*mse_fir(k)+(1-0.9735)*err_fir_e;
            mse_adf(k)=0.9735*mse_adf(k)+(1-0.9735)*err_adf_e;
            mse_mic(k)=0.9735*mse_mic(k)+(1-0.9735)*(abs(input_f(k)).^2);
            
            if mse_adf(k) > mse_mic(k) * 8
                adf_coeff_hi(k-32,:)=0;
                mse_adf(k)=mse_mic(k);
                err_adf=input_f(k);
            else
                if mse_adf(k) < min(0.125 * mse_mic(k), 0.5 * mse_fir(k)) && ref_cnt == 0
                    fir_coeff_hi(k-32,:)=adf_coeff_hi(k-32,:);
                    mse_fir(k)=mse_adf(k);
                    err_fir=err_adf;
                end
            end
            
            if mse_fir(k) > mse_mic(k) * 8 && ref_cnt == 0
                fir_coeff_hi(k-32,:)=0;
                mse_fir(k)=mse_mic(k);
                err_fir=input_f(k);
            else
                if mse_fir(k) < min(0.125 * mse_mic(k), 0.5 * mse_adf(k))
                    adf_coeff_hi(k-32,:)=fir_coeff_hi(k-32,:);
                    mse_adf(k)=mse_fir(k);
                    err_adf=err_fir;
                end
            end
            
            if k <= BAND_TABLE(6)
                subband_index=3;
            else
                subband_index=4;
            end
            
            if dt_flag == 2
                tmp = 80;
            else
                tmp = 20;
            end
            
            B=noise_est_mic_parameters(k).noise_level(1);
            if input_power * erl_ratio(subband_index) > tmp * B && filter_freeze == 0
                adf_coeff_hi(k-32,:) = adf_coeff_hi(k-32,:) + norm_step * err_adf' * stack_ref_hi(k-32,:);
                fir_update_flag(k)=1;
            else
                fir_update_flag(k)=0;
            end
            
            adf_early_tmp=input_f(k)-est_adf_early;
            fir_early_tmp=input_f(k)-est_fir_early;
            
            if abs(fir_early_tmp)^2 > abs(adf_early_tmp)^2
                min_power_early=abs(adf_early_tmp)^2;
                early_out=adf_early_tmp;
            else
                min_power_early=abs(fir_early_tmp)^2;
                early_out=fir_early_tmp;
            end
            
            if err_fir_e > err_adf_e
                power_ori_mic(k)=err_adf_e;
                power_echo(k)=est_adf_e;
                if err_adf_e > min_power_early
                    fir_out(k)=early_out;
                else
                    fir_out(k)=err_adf;
                end 
            else
                power_ori_mic(k)=err_fir_e;
                power_echo(k)=est_fir_e;
                if err_fir_e > min_power_early
                    fir_out(k)=early_out;
                else
                    fir_out(k)=err_fir;
                end
            end
            fir_est_ref(k)=input_f(k)-fir_out(k);         
        end
        for k=127:128
            fir_out(k)=0;
        end

        total_input_f(1:129)=fir_out(1:129);
        stack_est_ref(:,2:6) = stack_est_ref(:, 1:5);
        stack_est_ref(:, 1) = fir_est_ref.';
        
        %  erl
        input_mic_energy=zeros(1,4);
        est_mic_energy=zeros(1,4);
        echo_return_energy=zeros(1,4);
        input_f_e = abs(input_f).^2;
        for j=1:4
            input_mic_energy(j) = sum(input_f_e(BAND_TABLE(2*j-1):BAND_TABLE(2*j)));
            est_mic_energy(j) = sum(power_ori_mic(BAND_TABLE(2*j-1):BAND_TABLE(2*j)));
            echo_return_energy(j) = sum(power_echo(BAND_TABLE(2*j-1):BAND_TABLE(2*j)));
            
            if erl_peak(j) < echo_return_energy(j)
                erl_peak(j)=echo_return_energy(j);
            else
                erl_peak(j)=erl_peak(j)*0.9048+echo_return_energy(j)*(1-0.9048);
            end
            
            C = energy_group_parameters(j).noise_level(2);
            if energy_group(j) > 10 * C && input_mic_energy(j) > 4 * est_mic_energy(j) && filter_freeze == 0 && far_end_hold_time == 20
                erl_val = erl_peak(j) / (energy_group_peak(j) + 10^-6);
                erl_val = min(32, erl_val);
                erl_val = max(0.00001, erl_val);
                erl_ratio(j) = 0.8 * erl_ratio(j) + (1-0.8)*erl_val;
            end
        end
        
        min_erl_ratio = erl_ratio(2) / 8;
        max_erl_ratio = erl_ratio(2) * 8;
        
        for j=1:4
            erl_ratio(j)=max(erl_ratio(j), min_erl_ratio);
            erl_ratio(j)=min(erl_ratio(j), max_erl_ratio);
        end
        
        % nl_coeff1 is very much close to 1
        [nlp_parameters, nl_coeff1, fir_out] = nlp(nlp_parameters, ref_peak, fir_out, input_f_e, input_f, input_r, far_end_hold_time, WEB_RTC_AEC_NL_WEIGHT_CURVE);
        fir_out_e = abs(fir_out).^2;
        
        % double_talk
        % fir_out_e, ref_e, erl_ratio, B, far_end_hold_time
        dt_mic_psd_sum=0;
        mic_spk_ratio=0;
        for j=1:37
            freq=(j-1+0.5)*200 + 600;
            dt_spk_psd=sum(ref_e(DOUBLETALK_BAND_TABLE(2*j-1):DOUBLETALK_BAND_TABLE(2*j)-1));
            dt_mic_psd=sum(fir_out_e(DOUBLETALK_BAND_TABLE(2*j-1):DOUBLETALK_BAND_TABLE(2*j)-1));
            dt_ns_psd=0;
            for p=DOUBLETALK_BAND_TABLE(2*j-1):DOUBLETALK_BAND_TABLE(2*j)-1
                B=noise_est_mic_parameters(p).noise_level(1);
                dt_ns_psd = dt_ns_psd + B;
            end
            
            if freq > 800
                dt_mic_psd_sum = dt_mic_psd_sum + dt_mic_psd;
            end
            
            dt_mic_psd = dt_mic_psd - 10 * dt_ns_psd;
            dt_mic_psd = max(dt_mic_psd, 0);
            
            if dt_spk_peak(j) < dt_spk_psd
                dt_spk_peak(j) = dt_spk_psd;
            else
                dt_spk_peak(j) = 0.7 * dt_spk_peak(j) + 0.3 * dt_spk_psd;
            end
            
            if freq >= 0 && freq < 600 
                index = 1;
            elseif freq >= 600 && freq < 1200
                index = 2;
            elseif freq >= 1200 && freq < 3000
                index = 3;
            else
                index = 4;
            end
            
            if far_end_hold_time == 20
                mic_spk_ratio = mic_spk_ratio + dt_mic_psd/(dt_spk_peak(j) * erl_ratio(index)); % sum of ratio?
            end
        end
        
        if far_end_hold_time == 20
            [doubletalk_1_parameters] = aec_noise_estimation(dt_mic_psd_sum, doubletalk_1_parameters);
        end
        
        D=doubletalk_1_parameters.noise_level(1);
        mic_snr = dt_mic_psd_sum / D;
        if dt_frame_count < 10
            mic_ratio_smooth = mic_spk_ratio;
        else
            if mic_ratio_smooth < mic_spk_ratio
                mic_ratio_smooth = 0.95 * mic_ratio_smooth + 0.05 * mic_spk_ratio;
            else
                mic_ratio_smooth = 0.85 * mic_ratio_smooth + 0.15 * mic_spk_ratio;
            end
        end
        
        if far_end_hold_time == 20
            [doubletalk_2_parameters] = aec_noise_estimation(mic_ratio_smooth, doubletalk_2_parameters);
        end
        E=doubletalk_2_parameters.noise_level(1);
        mic_ratio_level = mic_spk_ratio / E;
        dt_flag = 0;
        if dt_frame_count < 10
            dt_frame_count = dt_frame_count +1;
        else
            if mic_ratio_level > 10 && mic_snr > 50
                dt_count =40;
            else
                dt_count = max(0, dt_count-1);
            end
            
            if far_end_hold_time == 0
                dt_flag = 1;
            else
                if dt_count > 0
                    dt_flag = 2;
                else
                    dt_flag = 0;
                end
            end
        end


        if dt_flag == 0
            for k = 3:125
                [noise_est_mic_parameters(k)] = aec_noise_estimation(power_ori_mic(k), noise_est_mic_parameters(k));
            end
        end
    else
        total_input_f(1:129) = input_f_bak;
    end

    % deal with 2nd mic, should not used if not for efficiency purpose
    input_f2=x2_f(i,:);
    input_f2_bak=input_f2;

    if no_ref_count < 500
        input_f2(1:4)=input_f2(1:4).*AEC_HPF_COEFF;
        %retf process
        % input_f2, stack_est_ref, fir_out2, fir_update_flag
        input_f2_e = abs(input_f2).^2;
        for k =1:129
            est_ref_fir = stack_est_ref(k,:) * fir_coeff_2(k, :)';
            est_ref_adf = stack_est_ref(k,:) * adf_coeff_2(k, :)';
            err_fir = input_f2(k) - est_ref_fir;
            err_adf = input_f2(k) - est_ref_adf;
            mu = 0.5 / (sum(abs(stack_est_ref(k, :).^2)) + 0.01);
            err_fir_e = abs(err_fir).^2;
            err_adf_e = abs(err_adf).^2;
            mse_fir2(k) = 0.9608 * mse_fir2(k) + (1-0.9608) * err_fir_e;
            mse_adf2(k) = 0.9608 * mse_adf2(k) + (1-0.9608) * err_adf_e;
            mse_mic2(k) = 0.9608 * mse_mic2(k) + (1-0.9608) * input_f2_e(k);
            
            if mse_adf2(k) > mse_mic2(k) * 8
                adf_coeff_2(k, :) = 0;
                mse_adf2(k) = mse_mic2(k);
                err_adf = input_f2(k);
            elseif mse_mic2(k) > mse_adf2(k) * 8 && mse_adf2(k) < 0.5 * mse_fir2(k)
                fir_coeff_2(k,:)= adf_coeff_2(k,:);
                mse_fir2(k) = mse_adf2(k);
                err_fir=err_adf;
            end
            
            if mse_fir2(k) > mse_mic2(k) * 8
                fir_coeff_2(k,:)=0;
                mse_fir2(k)=mse_mic2(k);
                err_fir =input_f2(k);
            elseif mse_mic2(k) > mse_fir2(k) * 8 && mse_fir2(k) < 0.5 * mse_adf2(k)
                adf_coeff_2(k,:) = fir_coeff_2(k,:);
                mse_adf2(k) = mse_fir2(k);
                err_adf=err_fir;
            end
            
            if fir_update_flag(k) == 1
                adf_coeff_2(k,:) = adf_coeff_2(k,:) + mu * err_adf' * stack_est_ref(k,:);
            end
            
            if err_adf_e <= err_fir_e && err_adf_e < input_f2_e(k)
                fir_out2(k) = err_adf; 
            elseif err_fir_e < err_adf_e && err_fir_e < input_f2_e(k)
                fir_out2(k) = err_fir;
            else
                fir_out2(k) = input_f2(k);
            end      
        end

        fir_out2 = fir_out2 .* min(1, abs(total_input_f(1:129)) .* abs(input_f2_bak) ./ (abs(fir_out2) .* abs(input_f_bak) + 10^-10));    
        total_input_f(130:258) = fir_out2;
        
        input_f_e2 = abs(input_f2_bak).^2;
        [nlp_parameters2, nl_coeff2, ~] = nlp(nlp_parameters2, ref_peak, fir_out2, input_f_e2, input_f2_bak, input_r, far_end_hold_time, WEB_RTC_AEC_NL_WEIGHT_CURVE); 
    else
        total_input_f(130:258)=input_f2_bak;
    end

    if no_ref_count < 500
        nl_coeff=mean([nl_coeff1;nl_coeff2]);
        nlp_overdrive = mean([nlp_parameters.over_drive, nlp_parameters2.over_drive]);
        
        % post_process
        % nlp_overdrive, total_input_f, nl_coeff, dt_flag, dt_st_strict
        nl_coeff_mean=mean(nl_coeff(5:24));
        [post_parameters] = aec_noise_estimation(nl_coeff_mean, post_parameters);
        G=post_parameters.noise_level(1);
        nlp_snr = nl_coeff_mean/(G + 10^-10);
        if dt_flag ~= 0 && nlp_snr > 3
            dt_hold_age = 40;
        else
            dt_hold_age = max(dt_hold_age-1, 0);
        end
        
        if dt_flag ~= 1
            if dt_hold_age > 0
                dt_flag = 2;
            else
                dt_flag = 0;
            end
        end
        
        over_drive_tmp = 0.25 * nlp_overdrive / (nlp_snr + 10^-10);
        post_over_drive_sm = 0.95 * post_over_drive_sm + 0.05 * over_drive_tmp;
        
        if dt_flag == 2
            over_drive_tmp = min(1, post_over_drive_sm);
        else
            over_drive_tmp = min(10, post_over_drive_sm);
        end

        tmp_ = nl_coeff .^ over_drive_tmp;
        if dt_flag == 2
           tmp_ = max(tmp_, 0.1);
        elseif dt_flag == 0
           tmp_ = max(tmp_, 0.01);
        end 
        
        total_input_f = total_input_f .* [tmp_, tmp_];  
    end
    
    if no_ref_count > 20
        dt_flag = 1;
    end

    total_input_f_pcm(i, :) = total_input_f;
 
end
close(hh);
end
