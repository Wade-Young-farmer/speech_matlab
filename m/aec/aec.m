function total_input_f_pcm=aec(mic_name, mic_name_2, ref_name, filter_coeff, WEB_RTC_AEC_NL_WEIGHT_CURVE, DOUBLETALK_BAND_TABLE, BAND_TABLE)
[r_enframe, r_f]=hpf(ref_name, filter_coeff);
[x_enframe, x_f]=hpf(mic_name, filter_coeff);
[x2_enframe, x2_f]=hpf(mic_name_2, filter_coeff);

ref_peak=0;
peak_hold_frame=0;

detect_far_end_buffer=zeros(1, 6);
% far_end_talk_flag=0;
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
dt_flag_strict=0;
fir_update_flag=zeros(1, 129);
fir_out=zeros(1, 129);
fir_out2=zeros(1, 129);
power_ori_mic=zeros(1, 129);
fir_est_ref=zeros(1, 129);
power_echo=zeros(1,129);
% used in fir part, strang meaning
ref_cnt=0;

erl_peak=zeros(1,4);
res_psd_smooth=zeros(1,129);
near_psd_smooth=zeros(1,129);
ref_psd_smooth=zeros(1,129);
corr_near_res_smooth=zeros(1,129);
corr_near_ref_smooth=zeros(1,129);
nlp_coeff_temp=1.0;
near_state_hold=0;
near_state=0;
over_drive=2;
over_drive_smooth=2;

nl_coeff_fb_min=0;
nl_coeff_fb_local_min=1;
nlp_is_new_min=0;
nlp_new_min_ctrl=0;
nlp_min_hold_time=0;
min_coeff_tmp=1;

res_psd_smooth2=zeros(1,129);
near_psd_smooth2=zeros(1,129);
ref_psd_smooth2=zeros(1,129);
corr_near_res_smooth2=zeros(1,129);
corr_near_ref_smooth2=zeros(1,129);
nlp_coeff_temp2=1.0;
near_state_hold2=0;
near_state2=0;
over_drive2=2;
over_drive_smooth2=2;

nl_coeff_fb_min2=0;
nl_coeff_fb_local_min2=1;
nlp_is_new_min2=0;
nlp_new_min_ctrl2=0;
nlp_min_hold_time2=0;
min_coeff_tmp2=1;

dt_spk_peak=zeros(1, 37);
dt_frame_count = 0;
mic_ratio_smooth = 0;
dt_count=0;
dt_count_strict=0;

total_input_f=zeros(1, 258);
stack_est_ref=zeros(129,6);
dt_hold_age=40;
post_over_drive_sm=0;

fir_coeff_2=zeros(129,6);
adf_coeff_2=zeros(129,6);
mse_fir2=zeros(1, 129);
mse_adf2=zeros(1, 129);
mse_mic2=zeros(1, 129);

spk_t_ns_parameters = initial(16, 100, 200);
for i=1:4
    energy_group_parameters(i)=initial(0.01, 16, 200);
end

for i=1:129
    noise_est_mic_parameters(i)=initial(0.1, 100, 200);
end
doubletalk_1_parameters=initial(0.1, 100, 62);
doubletalk_2_parameters=initial(0.01, 4, 200);
post_parameters=initial(0.01, 1, 100);

hh=waitbar(0, 'data is being processed');
pcm_size=size(x_enframe,1);
total_input_f_pcm=zeros(pcm_size, 258);
fir_out_debug=zeros(pcm_size, 10);
nl_coeff1_debug=zeros(pcm_size, 10);
dt_flag_debug=zeros(pcm_size,2);
total_input_1_f_debug=zeros(pcm_size, 10);
for i=1:size(x_enframe,1)
    str=[num2str(i), ' / ', num2str(pcm_size), ' processed']; 
    waitbar(i/pcm_size, hh, str); 
    
    input_t=x_enframe(i,:);
%     input_t2=x2_enframe(i,:);
    input_f=x_f(i,:);
    input_f2=x2_f(i,:);
    input_f_bak=input_f;
    input_f2_bak=input_f2;
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
    % [detect_far_end_buffer(1), spk_t_ns_parameters.noise_level]
    if detect_far_end_buffer(1) > max(10, 2*A)
%         far_end_talk_flag = 1;
        far_end_hold_time = 20;
%         i
    else
%         far_end_talk_flag = 0;
        far_end_hold_time = max(0, far_end_hold_time - 1);
    end
%     far_end_hold_time
    if far_end_hold_time == 0 && ref_peak < 20
        no_ref_count=min(1000, no_ref_count+1);
    else
        no_ref_count=0;
    end
%     no_ref_count
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
%     energy_group
    
    stack_ref_low(:, 2:20)=stack_ref_low(:, 1:19);
    stack_ref_low(:, 1)=ref_f(1:32).';
    stack_ref_hi(:,2:16)=stack_ref_hi(:,1:15);
    stack_ref_hi(:,1)=ref_f(33:128).';
    
%     no_ref_count
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
            if k == 4
                [err_fir_e, err_adf_e, min_power_early]
            end
            if err_fir_e > err_adf_e
                power_ori_mic(k)=err_adf_e;
                power_echo(k)=est_adf_e;
                if err_adf_e > min_power_early
                    fir_out(k)=early_out;
%                     if k == 4
%                         1
%                     end
                else
                    fir_out(k)=err_adf;
%                     if k == 4
%                         2
%                     end
                end 
            else
                power_ori_mic(k)=err_fir_e;
                power_echo(k)=est_fir_e;
                if err_fir_e > min_power_early
                    fir_out(k)=early_out;
%                     if k == 4
%                         3
%                     end
                else
                    fir_out(k)=err_fir;
%                     if k == 4
%                         4
%                     end
                end
            end
            fir_est_ref(k)=input_f(k)-fir_out(k);
            
%             if k == 4
%                 [early_out, err_adf, err_fir, fir_out(k)]                
%             end
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
        fir_out_debug(i,1:5)=fir_out(1:5);
        fir_out_debug(i,6:10) = fir_out(125:129);
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
              
%       erl_ratio
        
        % nlp
        % input_r, fir_out, input_f, ref_peak?????ar_end_hold_time, 1, 0
        if ref_peak > 5000
            volumn = 1;
        else
            volumn = 0;
        end
        
%         volumn
        res_psd=abs(fir_out).^2;
        near_psd=input_f_e;
        ref_psd=abs(input_r).^2;
        res_psd_smooth=0.7165*res_psd_smooth + (1-0.7165)*res_psd;
        near_psd_smooth=0.7165*near_psd_smooth + (1-0.7165)*near_psd;
        ref_psd_smooth=0.7165*ref_psd_smooth + (1-0.7165)*max(16,ref_psd);
        corr_near_res=conj(input_f) .* fir_out;
        corr_near_res_smooth=0.7165*corr_near_res_smooth + (1-0.7165)*corr_near_res;
        coh_temp_1=(abs(corr_near_res_smooth).^2) ./ (res_psd_smooth .* near_psd_smooth + 10^-10);
        corr_near_ref=conj(input_f) .* input_r;
        corr_near_ref_smooth=0.7165*corr_near_ref_smooth + (1-0.7165)*corr_near_ref;
        coh_temp_2=(abs(corr_near_ref_smooth).^2) ./ (near_psd_smooth .* ref_psd_smooth + 10^-10);
        
        coh_near_ref_avg=mean(coh_temp_2(5:32));
        coh_near_res_avg=mean(coh_temp_1(5:32)); % near is input here
        
        if 1 - coh_near_ref_avg < min(0.75, nlp_coeff_temp)
            nlp_coeff_temp = 1 - coh_near_ref_avg;
        end
        
        if coh_near_res_avg > 0.8 && coh_near_ref_avg < 0.3
            near_state = 1;           
            near_state_hold = 0;
        elseif coh_near_res_avg < 0.75 || coh_near_ref_avg > 0.5
            if near_state_hold == 3
                near_state = 0;
            else
                near_state_hold = near_state_hold + 1;
            end
        end
        
        if far_end_hold_time == 0
            near_state = 1; % represents someone is talking or speaker is not talking
            near_state_hold = 0;
        end
        
        if nlp_coeff_temp == 1 % represents no far end
            over_drive=1;
            if near_state == 1
                nl_coeff=coh_temp_1;
                nl_coeff_fb = coh_near_res_avg;
                nl_coeff_fb_low = coh_near_res_avg;
            else
                nl_coeff=1- coh_temp_2;
                nl_coeff_fb=1- coh_near_ref_avg;
                nl_coeff_fb_low= 1- coh_near_ref_avg;
            end
        else
            if volumn == 0 && near_state == 1
                nl_coeff=coh_temp_1;
                nl_coeff_fb = coh_near_res_avg;
                nl_coeff_fb_low = coh_near_res_avg;
            else
                nl_coeff=min(coh_temp_1, 1-coh_temp_2);
                nl_coeff_temp_array=sort(nl_coeff(5:32));
                nl_coeff_fb = nl_coeff_temp_array(21);
                nl_coeff_fb_low = nl_coeff_temp_array(14);
            end
        end
        
        nlp_coeff_temp = min(nlp_coeff_temp+0.0003, 1);
        
        % Min track nl_coeff_fb is not used in none wakeup judge case
        if nl_coeff_fb_low < min(0.6, nl_coeff_fb_local_min)
            nl_coeff_fb_min = nl_coeff_fb_low;
            nl_coeff_fb_local_min = nl_coeff_fb_low;
            nlp_is_new_min=1;
            nlp_new_min_ctrl=0;
            nlp_min_hold_time=0;
        else
            nlp_min_hold_time=nlp_min_hold_time+1;
        end
        
        if nlp_min_hold_time > 100 && nl_coeff_fb_low < min_coeff_tmp
            min_coeff_tmp = nl_coeff_fb_low;
        end
           
        if nlp_min_hold_time > 300 
            nl_coeff_fb_local_min = min_coeff_tmp;
            nl_coeff_fb_min = min_coeff_tmp;
            min_coeff_tmp = 1;
            nlp_min_hold_time = 150;
        end
        
        nl_coeff_fb_local_min = min(nl_coeff_fb_local_min + 0.0004, 1);
        if nlp_is_new_min == 1
            nlp_new_min_ctrl = nlp_new_min_ctrl + 1;
        end
        if nlp_new_min_ctrl == 2
            nlp_is_new_min = 0;
            nlp_new_min_ctrl = 0;
            over_drive = max(-1.15 / (log(nl_coeff_fb_min + 10^-10) + 10^-10), 1);
        end
        
        if over_drive < over_drive_smooth
            over_drive_smooth = 0.99 * over_drive_smooth + 0.01 * over_drive;
        else
            over_drive_smooth = 0.9 * over_drive_smooth + 0.1 * over_drive;
        end
  
%         nl_coeff = abs(nl_coeff); 
        for k=3:126
            if nl_coeff(k) > nl_coeff_fb
               nl_coeff(k) = (1-WEB_RTC_AEC_NL_WEIGHT_CURVE(k)) * nl_coeff(k) + WEB_RTC_AEC_NL_WEIGHT_CURVE(k) * nl_coeff_fb;
            end
            
            % judge difference is here! If there is false wakeup judge,
            % should multiply a different value
            fir_out(k) = fir_out(k) * nl_coeff(k);
        end
        fir_out_debug(i,1:5)=fir_out(1:5);
        fir_out_debug(i,6:10) = fir_out(125:129);
        % nl_coeff1 is very much close to 1
        nl_coeff1 = nl_coeff;
        nl_coeff1_debug(i,1:5)=nl_coeff1(1:5);
        nl_coeff1_debug(i,6:10)=nl_coeff1(125:129);
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
        dt_flag_strict = 0;
        if dt_frame_count < 10
            dt_frame_count = dt_frame_count +1;
        else
            if mic_ratio_level > 10 && mic_snr > 50
                dt_count =40;
            else
                dt_count = max(0, dt_count-1);
            end
            
            if dt_count > 0
                dt_flag = 2;
            end
            
            if far_end_hold_time == 0
                dt_flag = 1;
            end
            
            if mic_ratio_level > 50 && mic_snr > 100
                dt_count_strict = 40;
            else
                dt_count_strict = max(0, dt_count_strict -1);
            end
            
            if dt_count_strict > 0
                dt_flag_strict = 2;
            end
            
            if far_end_hold_time == 0
                dt_flag_strict = 1;
            end
        end
%         dt_frame_count
%         far_end_hold_time
        dt_flag_debug(i,:) = [dt_flag, dt_flag_strict];
        if dt_flag == 0
            for k = 3:125
                [noise_est_mic_parameters(k)] = aec_noise_estimation(power_ori_mic(k), noise_est_mic_parameters(k));
            end
        end
    else
        total_input_f(1:129) = input_f_bak;
    end
    
    total_input_1_f_debug(i, 1:5) = total_input_f(1:5);
    total_input_1_f_debug(i, 6:10) = total_input_f(125:129);
    
    % deal with 2nd mic, should not used if not for efficiency purpose    
    if no_ref_count < 500
        input_f2(1:4)=input_f2(1:4).*AEC_HPF_COEFF;
        %retf process
        % input_f2, stack_est_ref, fir_out2, fir_update_flag
        input_f2_e = abs(input_f2).^2;
        for k =1:129
            est_ref_fir = stack_est_ref(k,:) * fir_coeff_2(k, :)';
            est_ref_adf = stack_est_ref(k,:) * adf_coeff_2(k, :)';
%             est_ref_fir_e = abs(est_ref_fir).^2;
%             est_ref_adf_e = abs(est_ref_adf).^2;
            err_fir = input_f2(k) - est_ref_fir;
            err_adf = input_f2(k) - est_ref_adf;
            mu = 0.5 / (sum(abs(stack_est_ref(k, :).^2) + 0.01));
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
            
            if err_fir_e >= err_adf_e && err_adf_e < input_f2_e(k)
                fir_out2(k) = err_adf;
            elseif err_fir_e >= err_adf_e && err_adf_e >= input_f2_e(k)
                fir_out2(k) = input_f2(k); 
            elseif err_fir_e < err_adf_e && err_fir_e < input_f2_e(k)
                fir_out2(k) = err_fir;
            else
                fir_out2(k) = input_f2(k);
            end      
        end
        fir_out2 = fir_out2 .* max(1, abs(total_input_f(1:129)) .* abs(input_f2_bak) / (abs(fir_out2) .* abs(input_f_bak) + 10^-10));    
        total_input_f(130:258) = fir_out2;
        
        % nlp
        % use input_f2_bak instead of input_f2, should be the same
        % input_r, fir_out2, input_f2_bak, ref_peak?????ar_end_hold_time, 0, 0
        if ref_peak > 5000
            volumn = 1;
        else
            volumn = 0;
        end        
        res_psd=abs(fir_out2).^2;
        near_psd=abs(input_f2).^2;
        ref_psd=abs(input_r).^2;
        res_psd_smooth2=0.7165*res_psd_smooth2 + (1-0.7165)*res_psd;
        near_psd_smooth2=0.7165*near_psd_smooth2 + (1-0.7165)*near_psd;
        ref_psd_smooth2=0.7165*ref_psd_smooth2 + (1-0.7165)*max(16,ref_psd);
        corr_near_res=conj(input_f2) .* fir_out2;
        corr_near_res_smooth2=0.7165*corr_near_res_smooth2 + (1-0.7165)*corr_near_res;
        coh_temp_1=(abs(corr_near_res_smooth2).^2) ./ (res_psd_smooth .* near_psd_smooth + 10^-10);
        % res_psd_temp=sum(res_psd(3:126));
        % near_psd_temp=sum(near_psd(3:126));
        corr_near_ref=conj(input_f) .* input_r;
        corr_near_ref_smooth2=0.7165*corr_near_ref_smooth2 + (1-0.7165)*corr_near_ref;
        coh_temp_2=(abs(corr_near_ref_smooth2).^2) ./ (near_psd_smooth .* ref_psd_smooth + 10^-10);
        
        coh_near_ref_avg=mean(coh_temp_2(5:32));
        coh_near_res_avg=mean(coh_temp_1(5:32)); % near is input here
        
        if 1 - coh_near_ref_avg < min(0.75, nlp_coeff_temp2)
            nlp_coeff_temp2 = 1 - coh_near_ref_avg;
        end
        
        if coh_near_res_avg > 0.8 && coh_near_ref_avg < 0.3
            near_state2 = 1; % ?????????state????????????????????????????????????????lp????????????            
            near_state_hold2 = 0;
        elseif coh_near_res_avg < 0.75 || coh_near_ref_avg > 0.5
            if near_state_hold2 == 3
                near_state2 = 0;
            else
                near_state_hold2 = near_state_hold2 + 1;
            end
        end
        
        if far_end_hold_time == 0
            near_state2 = 1; % represents someone is talking or speaker is not talking
            near_state_hold2 = 0;
        end
        
        if nlp_coeff_temp2 == 1 % represents no far end
            over_drive2=1;
            if near_state2 == 1
                nl_coeff=coh_temp_1;
                nl_coeff_fb = coh_near_res_avg;
                nl_coeff_fb_low = coh_near_res_avg;
                % ??ear_res
            else
                % ??- near_far
                nl_coeff=1- coh_temp_2;
                nl_coeff_fb=1- coh_near_ref_avg;
                nl_coeff_fb_low= 1- coh_near_ref_avg;
            end
        else
            if volumn == 0 && near_state2 == 1
                nl_coeff=coh_temp_1;
                nl_coeff_fb = coh_near_res_avg;
                nl_coeff_fb_low = coh_near_res_avg;
                % ??near_res
            else
                % ??min(near_res, 1 - near_far)
                nl_coeff=min(coh_temp_1, 1-coh_temp_2);
                nl_coeff_temp_array=sort(nl_coeff(5:32));
                nl_coeff_fb = nl_coeff_temp_array(21);
                nl_coeff_fb_low = nl_coeff_temp_array(14);
            end
        end
        
        nlp_coeff_temp2 = min(nlp_coeff_temp2+0.0003, 1);
        
        % Min track nl_coeff_fb?????is not used in none wakeup judge case
        if nl_coeff_fb_low < min(0.6, nl_coeff_fb_local_min2)
            nl_coeff_fb_min2 = nl_coeff_fb_low;
            nl_coeff_fb_local_min2 = nl_coeff_fb_low;
            nlp_is_new_min2=1;
            nlp_new_min_ctrl2=0;
            nlp_min_hold_time2=0;
        else
            nlp_min_hold_time2=nlp_min_hold_time2+1;
        end
        
        if nlp_min_hold_time2 > 100 && nl_coeff_fb_low < min_coeff_tmp2
            min_coeff_tmp2 = nl_coeff_fb_low;
        end
           
        if nlp_min_hold_time2 > 300 
            nl_coeff_fb_local_min2 = min_coeff_tmp2;
            nl_coeff_fb_min2 = min_coeff_tmp2;
            min_coeff_tmp2 = 1;
            nlp_min_hold_time2 = 150;
        end
        
        nl_coeff_fb_local_min2 = min(nl_coeff_fb_local_min2 + 0.0004, 1);
        if nlp_is_new_min2 == 1
            nlp_new_min_ctrl2 = nlp_new_min_ctrl2 + 1;
        end
        if nlp_new_min_ctrl2 == 2
            nlp_is_new_min2 = 0;
            nlp_new_min_ctrl2 = 0;
            over_drive2 = max(-1.15 / (log(nl_coeff_fb_min2 + 10^-10) + 10^-10), 1);
        end
        
        if over_drive2 < over_drive_smooth2
            over_drive_smooth2 = 0.99 * over_drive_smooth2 + 0.01 * over_drive2;
        else
            over_drive_smooth2 = 0.9 * over_drive_smooth2 + 0.1 * over_drive2;
        end
        
%         nl_coeff = abs(nl_coeff); 
        for k=3:126
            if nl_coeff(k) > nl_coeff_fb
               nl_coeff(k) = (1 -WEB_RTC_AEC_NL_WEIGHT_CURVE(k)) * nl_coeff(k) + WEB_RTC_AEC_NL_WEIGHT_CURVE(k) * nl_coeff_fb;
            end
        end
        nl_coeff2 = nl_coeff;
    else
        total_input_f(130:258)=input_f2_bak;
    end
    
    if no_ref_count < 500
        nl_coeff=mean([nl_coeff1;nl_coeff2]);
        nlp_overdrive = mean([over_drive, over_drive2]);
        
        % post_process
        % nlp_overdrive, total_input_f, nl_coeff, dt_flag, dt_st_strict, 2
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
        
        for k = 1:129
            tmp_ = nl_coeff(k).^over_drive_tmp;
            if dt_flag == 2
                tmp_ = max(tmp_, 0.1);
            elseif dt_flag == 0
                tmp_ = max(tmp_, 0.01);
            end
            
            total_input_f(k) = total_input_f(k) * tmp_;
            total_input_f(k+129) = total_input_f(k+129) * tmp_;
        end  
    end
    
    if no_ref_count > 20
        dt_flag = 1;
    end
    total_input_f_pcm(i, :) = total_input_f;
        
    if i == 481
       save fir_out_debug.mat fir_out_debug
       save nl_coeff1_debug.mat nl_coeff1_debug
       save dt_flag_debug.mat dt_flag_debug
       save total_input_1_f_debug.mat total_input_1_f_debug
       pause;
    end
end
close(hh);
end
