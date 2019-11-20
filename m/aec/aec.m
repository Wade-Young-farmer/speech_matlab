function my_aec(mic_name, ref_name)
Srate=16000;
file_id=fopen(mic_name, 'r');
x=fread(file_id, inf, 'int16');
fclose(file_id);

file_id=fopen(ref_name, 'r');
r=fread(file_id, inf, 'int16');
fclose(file_id);

x_enframe = enframe(x,128);
input_buffer=zeros(1, 256);

r_enframe = enframe(r,128);
ref_buffer=zeros(1, 256);

ref_peak=0;
smooth_factor1=1-exp(-1/10);
peak_hold_frame=0;

detect_far_end_buffer=zeros(1, 6);
far_end_talk_flag=0;
far_end_hold_time=1;
no_ref_count=0;

BAND_TABLE= [2, 8, 9, 18, 19, 47, 48, 124];
BAND_TABLE=BAND_TABLE+1;
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
power_ori_mic=zeros(1, 129);
fir_est_ref=zeros(1, 129);
power_echo=zeros(1,129);

erl_peak=zeros(1,4);

for i=1:size(x_enframe,1)
    input_buffer(1:128)=input_buffer(129:256);
    input_buffer(129:256)=x_enframe(i,:);
    input_f=fft(input_buffer);
    input_f=input_f(1:129);
    input_t=x_enframe(i,:);
    
    ref_buffer(1:128)=ref_buffer(129:256);
    ref_buffer(129:256)=r_enframe(i,:);
    ref_f=fft(ref_buffer);
    ref_f=ref_f(1:129);
    
    % Only apply to single-channel reference signal
    ref_t=r_enframe(i,:);
    if max(abs(ref_t)) > ref_peak
        ref_peak = max(abs(ref_t));
        peak_hold_frame=0;
    else
        peak_hold_frame=peak_hold_frame+1;
        if peak_hold_frame > 100
            ref_peak=(1-smooth_factor1)*ref_peak + smooth_factor1*max(abs(ref_t));
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
    
    A=1; % A是一级ns的值, 由detect_far_end_buffer(1)得到
    if detect_far_end_buffer(1) > max(10, 2*A)
        far_end_talk_flag = 1;
        far_end_hold_time = 20;
    else
        far_end_talk_flag = 0;
        far_end_hold_time = min(0, far_end_hold_time - 1);
    end
    
    if far_end_hold_time == 0 && ref_peak < 20
        no_ref_count=max(1000, no_ref_count+1);
    else
        no_ref_count=0;
    end
        
    ref_e = abs(ref_f).^2;
    for j=1:4
        energy_group(j)=0;
        for k=BAND_TABLE(j):BAND_TABLE(j+1)
            energy_group(j)=energy_group(j)+ref_e(k);
        end
        if energy_group(j) > energy_group_peak(j)
            energy_group_peak(j) = energy_group(j);
        else
            energy_group_peak(j) = 0.9048 * energy_group_peak(j) + (1-0.9048) * energy_group(j);
        end
        
        C=1; % energy_group(j) 跟踪出一个ref子带域二级ns的值
    end
    
    stack_ref_low(:, 2:20)=stack_ref_low(:, 1:19);
    stack_ref_low(:, 1)=ref_f(1:32).';
    stack_ref_hi(:,2:16)=stack_ref_hi(:,1:15);
    stack_ref_hi(:,1)=ref_f(33:128).';
    
    if no_ref_count < 500
        input_f(1:4)=input_f(1:4).*AEC_HPF_COEFF;
        mic_peak=max(abs(input_t));
        if far_end_hold_time < 20
            ref_cnt=20; % naming not reasonable
        else
            ref_cnt=max(0, ref_cnt-1);
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
            est_fir=stack_ref_low(k,:) * fir_coeff_low(k,:)'; % Matrix Mul, waring of not a scalar
            est_adf=stack_ref_low(k,:) * adf_coeff_low(k,:)';
            est_fir_e = abs(est_fir).^2;
            est_adf_e = abs(est_adf).^2;
            err_fir=input_f(k)-est_fir;
            err_adf=input_f(k)-est_adf;
            input_power=sum(abs(stack_ref_low(k,:)).^2);
            norm_step=0.5/(input_power+0.01);
            err_fir_e=abs(err_fir).^2;
            err_adf_e=abs(err_adf).^2;
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
            
            B=1; % B是一级ns的值，每个mic的每个频点都会估计出一个
            if input_power * erl_ratio(subband_index) > tmp * B && filter_freeze == 0
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
                early_out=fir_early_out;
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
            est_fir_early=stack_ref_hi(k,1:4) * fir_coeff_hi(k,1:4)';
            est_adf_early=stack_ref_hi(k,1:4) * adf_coeff_hi(k,1:4)';
            est_fir=stack_ref_hi(k,:) * fir_coeff_hi(k,:)'; % Matrix Mul, waring of not a scalar
            est_adf=stack_ref_hi(k,:) * adf_coeff_hi(k,:)';
            est_fir_e = abs(est_fir).^2;
            est_adf_e = abs(est_adf).^2;
            err_fir=input_f(k)-est_fir;
            err_adf=input_f(k)-est_adf;
            input_power=sum(abs(stack_ref_hi(k,:)).^2);
            norm_step=0.5/(input_power+0.01);
            err_fir_e=abs(err_fir).^2;
            err_adf_e=abs(err_adf).^2;
            mse_fir(k)=0.9735*mse_fir(k)+(1-0.9735)*err_fir_e;
            mse_adf(k)=0.9735*mse_adf(k)+(1-0.9735)*err_adf_e;
            mse_mic(k)=0.9735*mse_mic(k)+(1-0.9735)*(abs(input_f(k)).^2);
            
            if mse_adf(k) > mse_mic(k) * 8
                adf_coeff_hi(k,:)=0;
                mse_adf(k)=mse_mic(k);
                err_adf=input_f(k);
            else
                if mse_adf(k) < min(0.125 * mse_mic(k), 0.5 * mse_fir(k)) && ref_cnt == 0
                    fir_coeff_hi(k,:)=adf_coeff_hi(k,:);
                    mse_fir(k)=mse_adf(k);
                    err_fir=err_adf;
                end
            end
            
            if mse_fir(k) > mse_mic(k) * 8 && ref_cnt == 0
                fir_coeff_hi(k,:)=0;
                mse_fir(k)=mse_mic(k);
                err_fir=input_f(k);
            else
                if mse_fir(k) < min(0.125 * mse_mic(k), 0.5 * mse_adf(k))
                    adf_coeff_hi(k,:)=fir_coeff_hi(k,:);
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
            
            B=1; % B是一级ns的值，每个mic的每个频点都会估计出一个
            if input_power * erl_ratio(subband_index) > tmp * B && filter_freeze == 0
                adf_coeff_hi(k,:) = adf_coeff_hi(k,:) + norm_step * err_adf' * stack_ref_hi(k,:);
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
                early_out=fir_early_out;
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
        
        %  erl
        input_mic_energy=zeros(1,4);
        est_mic_energy=zeros(1,4);
        echo_return_energy=zeros(1,4);
        for j=1:4
            for k=BAND_TABLE(j):BAND_TABLE(j+1)
                input_mic_energy(j)=input_mic_energy(j)+abs(input_f(k))^2;
                est_mic_energy(j)=est_mic_energy(j)+power_ori_mic(k);
                echo_return_energy(j)=echo_return_energy(j)+power_echo(k);
            end
        end
        
        for j=1:4
            if erl_peak(j)<echo_return_energy(j)
                erl_peak(j)=echo_return_energy(j);
            else
                erl_peak(j)=erl_peak(j)*0.9048+echo_return_energy(j)*(1-0.9048);
            end
            
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
        
        % nlp
        % input_r, fir_out, input_f, ref_peak
        
        
        
    end
    
    
    
    
    
    
        
    
end
end