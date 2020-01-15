% len_val = 512;
% test = random('unif', 0, 1, [1, 512]);
% ns_ps = fft(test).';
% 
% parameters = struct('n',2,'len',len_val,'ad',0.95,'as',0.8,'ap',0.2,'beta',0.8,'beta1',0.98,'gamma',0.998, 'delta',2, 'alpha',0.7,...
%             'pk',zeros(len_val,1),'noise_ps',ns_ps,'pxk_old',ns_ps,'pxk',ns_ps,'pnk_old',ns_ps,'pnk',ns_ps);
% 
% ss = mcra2_estimation(ns_ps,parameters);


NoiseMcra.first  = 1;
NoiseMcra.len    = 257;
% NoiseMcra.w    = 1;
% NoiseMcra.L    = 60;
NoiseMcra.alpha    = 0.9;
% NoiseMcra.frm_cnt    = 0;
NoiseMcra.S   = zeros(1, NoiseMcra.len);
NoiseMcra.Smin   = zeros(1, NoiseMcra.len);
% NoiseMcra.Stmp   = zeros(1, NoiseMcra.len);
NoiseMcra.Yprob   = zeros(1, NoiseMcra.len);
NoiseMcra.lamda_d   = zeros(1, NoiseMcra.len);
NoiseMcra.P = zeros(1, NoiseMcra.len);
% NoiseMcra.b   = zeros(1, 2*NoiseMcra.w+1);    %hanning window coefficients

test = random('unif', 0, 1, [1, 512]);
ns_ps = fft(test);

ns_ps = ns_ps(1:NoiseMcra.len);
ns_ps = (abs(ns_ps)).^2;

for i = 1: 10
    [noise, NoiseMcra, Sr] = mcra(ns_ps, NoiseMcra);
end
