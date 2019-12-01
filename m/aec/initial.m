function parameters=initial(min_val, max_val, win_len)
min_noise_energy = min_val;
max_noise_energy = max_val;

min_energy=ones(1, 2)*max_val;
temp_min=ones(1,2)* max_val;
min_hold_frame=zeros(1,2);
noise_level=ones(1,2) * max_val;
min_win_len=[2,1]*win_len;

parameters = struct('min_noise_energy', min_noise_energy, 'max_noise_energy', max_noise_energy, 'min_energy', min_energy, ...
                    'temp_min', temp_min, 'min_hold_frame', min_hold_frame, 'noise_level', noise_level, ...
                    'min_win_len', min_win_len);
end