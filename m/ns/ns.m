% This is for speech enhancement 
function [parameters, noise_psd] = ns(parameters, mic_psd)
frame_count = parameters.frame_count;
hist_buf = parameters.hist_buf;
p_run_min_power_smoothed = parameters.p_run_min_power_smoothed;
p_run_min_power = parameters.p_run_min_power;
p_run_min_tmp = parameters.p_run_min_tmp;
p_run_min_age = parameters.p_run_min_age;

frame_count = frame_count + 1;
mean(mic_psd);
if mean(mic_psd) < 128 * 1e-7
    noise_psd = p_run_min_power_smoothed * 4;
    noise_psd = max(noise_psd, 3.16e-7 * 128);
else
    hist_buf(2:7,:) = hist_buf(1:6,:);
    hist_buf(1, :) = mic_psd;
    if frame_count < 7
        p_run_min_power = mic_psd;
        p_run_min_power_smoothed = 0.75 * p_run_min_power_smoothed + 0.25 * p_run_min_power;
        noise_psd = p_run_min_power_smoothed * 4;
        noise_psd = max(noise_psd, 3.16e-7 * 128);
    else
        a = mean(hist_buf, 1);
        b = mean(a);
        
        if b < 3.16e-7 * 128
            noise_psd = p_run_min_power_smoothed * 4;
            noise_psd = max(noise_psd, 3.16e-7 * 128);
        else
            for i=1:129
                if p_run_min_power(i) > a(i)
                    p_run_min_power(i) = a(i);
                    p_run_min_tmp(i) = 1e10;
                    p_run_min_age(i) = 0;
                else
                    p_run_min_age(i) = p_run_min_age(i) + 1;
                end
                
                if p_run_min_age(i) >= 250 && p_run_min_tmp(i) > a(i)
                    p_run_min_tmp(i) = a(i);
                end
                
                if p_run_min_age(i) >= 500
                    p_run_min_power(i) = p_run_min_tmp(i);
                    p_run_min_tmp(i) = 1e10;
                    p_run_min_age(i) = 250;
                end
            end
            
            p_run_min_power_smoothed = 0.75 * p_run_min_power_smoothed + 0.25 * p_run_min_power;
            noise_psd = p_run_min_power_smoothed * 4;
            noise_psd = max(noise_psd, 3.16e-7 * 128); 
        end
    end
end

parameters.p_run_min_power_smoothed = p_run_min_power_smoothed;
parameters.hist_buf = hist_buf;
parameters.p_run_min_power = p_run_min_power;
parameters.frame_count = frame_count;
parameters.p_run_min_tmp = p_run_min_tmp;
parameters.p_run_min_age = p_run_min_age;
end