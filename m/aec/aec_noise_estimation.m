function vad = aec_noise_estimation(input, paremeters)
min_noise_energy = paremeters.min_noise_energy;
max_noise_energy = paremeters.max_noise_energy;

min_energy=paremeters.min_energy;
temp_min=paremeters.temp_min;
min_hold_frame=paremeters.min_hold_frame;
noise_level=paremeters.noise_level;
min_win_len=paremeters.min_win_len;

if min_energy(2) > input
    min_energy(2) = input;
    min_hold_frame(2) = 0;
    temp_min(2) = max_noise_energy;
else
    min_hold_frame(2) = min_hold_frame(2) + 1;
end

if min_hold_frame(2) > 0.5 * min_win_len(2) && temp_min(2) > input
    temp_min(2) = input;
end

if min_hold_frame(2) > 1.5 * min_win_len(2)
    min_energy(2) = temp_min(2);
    temp_min(2) = max_noise_energy;
    min_hold_frame(2) = 0.5 * min_win_len(2);
end

noise_level(2) = 0.8 * noise_level(2) + 0.2 * min_energy(2); 

if input < 10 * noise_level(2)
    input = max(input, min_noise_energy);
    if input < min_energy(1)
        min_energy(1) = input;
        min_hold_frame(1) = 0;
        temp_min(1) = max_noise_energy;
    else
        min_hold_frame(1) = min_hold_frame(1) + 1;
    end
    
    if min_hold_frame(1) > 0.5 * min_win_len(1) && temp_min(1) > input
        temp_min(1) = input;
    end
    
    if min_hold_frame(1) > 1.5 * min_win_len(1)
        min_energy(1) = temp_min(1);
        temp_min(1) = max_noise_energy;
        min_hold_frame(1) = 0.5 * min_win_len(1);
    end
    
    noise_level(1) = 0.8 * noise_level(1) + 0.2 * min_energy(1); 
end

if input > 20 * max(noise_level)
    vad = 1;
else
    vad = 0;
end

paremeters.min_energy = min_energy;
paremeters.temp_min = temp_min;
paremeters.min_hold_frame = min_hold_frame;
paremeters.noise_level = noise_level;
paremeters.min_win_len = min_win_len;
end