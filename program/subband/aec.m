function val=my_smooth(y,x,f)
val=f*x+(1-f)*y;

function [mic_out,mic_complex_out, dt_st, aec_status]=aec(aec_status,mic_complex_in,ref_complex_in,mic_in,ref_in,dt_st)
mic_out=0;
mic_complex_out=0;
dt_st=0;
aec_status=0;
%% TDE_ENABLE


frame_num = size(ref_in, 2);
for i = 1:frame_num
    %% track peak value of reference data in time-domain
    mr=max(abs(ref_in(:,i)));
    if aec_status.ref_peak_track < mr
        aec_status.ref_peak_track = mr;
        aec_status.peak_hold_frame=0;
    else
        aec_status.peak_hold_frame=aec_status.peak_hold_frame+1;
        if aec_status.peak_hold_frame > 100 % 最大值保持已经超过100帧, 开始平滑
            aec_status.ref_peak_track = my_smooth(aec_status.ref_peak_track,mr,peak_smooth_factor);
        end
    end
    
    %% detect the exsitence of reference signal
    % st_noise_est_spk_t
    if i == 1
    else
    end
    
    if aec_status.far_end_talk_flag == 1
        aec_status.far_end_holdtime = 20;
    else
        if aec_status.far_end_holdtime > 0
            aec_status.far_end_holdtime=aec_status.far_end_holdtime-1;
        end
    end
    
    if aec_status.far_end_holdtime == 0 && aec_status.ref_peak_track < 20 % 需要far_end_talk_flag为零已经持续20帧，并且最大值的限制
        if aec_status.no_ref_count < 1000
            aec_status.no_ref_count=aec_status.no_ref_count+1; %这个值最大取1000帧也就够了，再多也没有意义
        end
    else
        aec_status.no_ref_count=0;
    end
    
    %% 参考信号预处理
    ref_in=ref_complex_in(:,i);
    ref_in_energy=abs(ref_in).^2;
    ref_part_energy(1)=sum(ref_in_energy(3:9));
    ref_part_energy(2)=sum(ref_in_energy(10:19));
    ref_part_energy(3)=sum(ref_in_energy(20:48));
    ref_part_energy(4)=sum(ref_in_energy(49:125));
    
    for j=1:4
        if aec_status.ref_peak_energy(j) < ref_part_energy(j)
            aec_status.ref_peak_energy(j)=ref_part_energy(j);
        else
            aec_status.ref_peak_energy(j)=my_smooth(aec_status.ref_peak_energy(j),ref_part_energy(j),alpha_peak);
        end
        %% 估计频域的参考信号的噪声水平 st_noise_est_spk_subband
    end
    
    
    %% 参考信号开始处理
    % 分低频段和高频段， 每个频段有过往20帧该频点的信息
    aec_status.stack_sig_low(:,2:20)=aec_status.stack_sig_low(:,1:19);
    aec_status.stack_sig_low(:,1)=ref_in(1:32,i);
    aec_status.stack_sig_hi(:,2:16)=aec_status.stack_sig_hi(:,2:16);
    aec_status.stack_sig_hi(:,1)=ref_in(33:128,i);
    
    %% 麦克风信号开始处理
    mic_in_=mic_in(:,i);
    mic_complex_in_=mic_complex_in(:,i);
    mic_complex_in__=mic_complex_in_;
    
    if aec_status.no_ref_cnt < 500 % 不是一直无ref的情况，意思就是有可能有内噪，500 *0.008，就是说有内噪，或者播放结束4s内
        % 高频滤波
        mic_complex_in_(1:4)=mic_complex_in_(1:4).*high_pass_freq_coeff.';
        aec_status.fir_far_end_holdtime=aec_status.far_end_holdtime;
        
        % 线性fir
        % 输入 hand
        
    end
    
    
end
function [fir_output,fir_power_output]=aec_fir(aec_status,mic_in_,mic_complex_in_)
%aec_status包含stack_sig_low, stack_sig_hi, double_talk_result[0],noise_est_mic_chan
mm=max(abs(mic_in_));








