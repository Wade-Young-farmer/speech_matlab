clc;
clear;

NoiseMcra1.first  = 1;
NoiseMcra1.len    = 65;
NoiseMcra1.w    = 1;
NoiseMcra1.L    = 60;
NoiseMcra1.alpha    = 0.85;
NoiseMcra1.frm_cnt    = 0;
NoiseMcra1.S   = zeros(1, NoiseMcra1.len);
NoiseMcra1.Smin   = zeros(1, NoiseMcra1.len);
NoiseMcra1.Stmp   = zeros(1, NoiseMcra1.len);
NoiseMcra1.Yprob   = zeros(1, NoiseMcra1.len);
NoiseMcra1.lamda_d   = zeros(1, NoiseMcra1.len);
NoiseMcra1.b   = zeros(1, 2*NoiseMcra1.w+1);    %hanning window coefficients
NoiseMcra1.b = 0.5 * (1 - cos(2 * pi * (1:2*NoiseMcra1.w+1) / (2*NoiseMcra1.w+2)));
NoiseMcra1.b = NoiseMcra1.b ./ sum(NoiseMcra1.b);

aa = ones(1,65);
[noi_spec, NoiseMcra1] = mcra(aa,NoiseMcra1);