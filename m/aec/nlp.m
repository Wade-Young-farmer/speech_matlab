function [parameters, nl_coeff, fir_out]=nlp(parameters, ref_peak, fir_out, input_f_e, input_f, input_r, far_end_hold_time, WEB_RTC_AEC_NL_WEIGHT_CURVE)
res_psd_smooth=parameters.res_psd_smooth;
near_psd_smooth=parameters.near_psd_smooth;
ref_psd_smooth=parameters.ref_psd_smooth;
corr_near_res_smooth=parameters.corr_near_res_smooth;
corr_near_ref_smooth=parameters.corr_near_ref_smooth;
nlp_coeff_temp=parameters.nlp_coeff_temp;
near_state_hold=parameters.near_state_hold;
near_state=parameters.near_state;
over_drive=parameters.over_drive;
over_drive_smooth=parameters.over_drive_smooth;

nl_coeff_fb_min=parameters.nl_coeff_fb_min;
nl_coeff_fb_local_min=parameters.nl_coeff_fb_local_min;
nlp_is_new_min=parameters.nlp_is_new_min;
nlp_new_min_ctrl=parameters.nlp_new_min_ctrl;
nlp_min_hold_time=parameters.nlp_min_hold_time;
min_coeff_tmp=parameters.min_coeff_tmp;

if ref_peak > 5000
    volumn = 1;
else
    volumn = 0;
end
        
res_psd=abs(fir_out).^2;
ref_psd=abs(input_r).^2;
near_psd=input_f_e;

res_psd_smooth=0.7165*res_psd_smooth + (1-0.7165)*res_psd;
ref_psd_smooth=0.7165*ref_psd_smooth + (1-0.7165)*max(16,ref_psd);
near_psd_smooth=0.7165*near_psd_smooth + (1-0.7165)*near_psd;
        
corr_near_res=conj(input_f) .* fir_out;        
corr_near_res_smooth=0.7165*corr_near_res_smooth + (1-0.7165)*corr_near_res;        
coh_temp_1=(abs(corr_near_res_smooth).^2) ./ (res_psd_smooth .* near_psd_smooth + 10^-10);
        
corr_near_ref=conj(input_f) .* input_r;       
corr_near_ref_smooth=0.7165*corr_near_ref_smooth + (1-0.7165)*corr_near_ref;        
coh_temp_2=(abs(corr_near_ref_smooth).^2) ./ (ref_psd_smooth .* near_psd_smooth + 10^-10);
                
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
      nl_coeff=1 - coh_temp_2;
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
        
%         % Min track nl_coeff_fb is not used in none wakeup judge case
%         if nl_coeff_fb_low < min(0.6, nl_coeff_fb_local_min)
%             nl_coeff_fb_min = nl_coeff_fb_low;
%             nl_coeff_fb_local_min = nl_coeff_fb_low;
%             nlp_is_new_min=1;
%             nlp_new_min_ctrl=0;
%             nlp_min_hold_time=0;
%         else
%             nlp_min_hold_time=nlp_min_hold_time+1;
%         end
%         
%         if nlp_min_hold_time > 100 && nl_coeff_fb_low < min_coeff_tmp
%             min_coeff_tmp = nl_coeff_fb_low;
%         end
%            
%         if nlp_min_hold_time > 300 
%             nl_coeff_fb_local_min = min_coeff_tmp;
%             nl_coeff_fb_min = min_coeff_tmp;
%             min_coeff_tmp = 1;
%             nlp_min_hold_time = 150;
%         end
%         
%         nl_coeff_fb_local_min = min(nl_coeff_fb_local_min + 0.0004, 1);
%         if nlp_is_new_min == 1
%             nlp_new_min_ctrl = nlp_new_min_ctrl + 1;
%         end
%         if nlp_new_min_ctrl == 2
%             nlp_is_new_min = 0;
%             nlp_new_min_ctrl = 0;
%             over_drive = max(-1.15 / (log(nl_coeff_fb_min + 10^-10) + 10^-10), 1);
%         end
%         
%         if over_drive < over_drive_smooth
%             over_drive_smooth = 0.99 * over_drive_smooth + 0.01 * over_drive;
%         else
%             over_drive_smooth = 0.9 * over_drive_smooth + 0.1 * over_drive;
%         end
  
%         nl_coeff = abs(nl_coeff); 
for k=3:126
    if nl_coeff(k) > nl_coeff_fb
       nl_coeff(k) = (1 - WEB_RTC_AEC_NL_WEIGHT_CURVE(k)) * nl_coeff(k) + WEB_RTC_AEC_NL_WEIGHT_CURVE(k) * nl_coeff_fb;
    end
            
    % judge difference is here! If there is false wakeup judge,
    % should multiply a different value
    fir_out(k) = fir_out(k) * nl_coeff(k);
end

parameters.res_psd_smooth = res_psd_smooth;
parameters.near_psd_smooth = near_psd_smooth;
parameters.ref_psd_smooth = ref_psd_smooth;
parameters.corr_near_res_smooth = corr_near_res_smooth;
parameters.corr_near_ref_smooth = corr_near_ref_smooth;

parameters.nlp_coeff_temp = nlp_coeff_temp;
parameters.near_state_hold = near_state_hold;
parameters.near_state = near_state;
parameters.over_drive = over_drive;
parameters.over_drive_smooth = over_drive_smooth;

parameters.nl_coeff_fb_min = nl_coeff_fb_min;
parameters.nl_coeff_fb_local_min = nl_coeff_fb_local_min;
parameters.nlp_is_new_min = nlp_is_new_min;
parameters.nlp_new_min_ctrl = nlp_new_min_ctrl;
parameters.nlp_min_hold_time = nlp_min_hold_time;
parameters.min_coeff_tmp = min_coeff_tmp;
end