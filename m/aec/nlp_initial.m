function parameters=nlp_initial()
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

parameters = struct('res_psd_smooth', res_psd_smooth, 'near_psd_smooth', near_psd_smooth, 'ref_psd_smooth', ref_psd_smooth, ...
                    'corr_near_res_smooth', corr_near_res_smooth, 'corr_near_ref_smooth', corr_near_ref_smooth, ...
                    'nlp_coeff_temp', nlp_coeff_temp, 'near_state_hold', near_state_hold, ...
                    'near_state', near_state, 'over_drive', over_drive, 'over_drive_smooth', over_drive_smooth, ...
                    'nl_coeff_fb_min', nl_coeff_fb_min, 'nl_coeff_fb_local_min', nl_coeff_fb_local_min, ...
                    'nlp_is_new_min', nlp_is_new_min, 'nlp_new_min_ctrl', nlp_new_min_ctrl, ...
                    'nlp_min_hold_time', nlp_min_hold_time, 'min_coeff_tmp', min_coeff_tmp);
end