function [noise, NoiseMcra] = noise_mcra(energy, NoiseMcra)

% NoiseMcra.first  = 1;
% NoiseMcra.len    = 129;
% NoiseMcra.w    = 1;
% NoiseMcra.L    = 60;
% NoiseMcra.alpha    = 0.9;
% NoiseMcra.frm_cnt    = 0;
% NoiseMcra.S   = zeros(1, NoiseMcra.len);
% NoiseMcra.Smin   = zeros(1, NoiseMcra.len);
% NoiseMcra.Stmp   = zeros(1, NoiseMcra.len);
% NoiseMcra.Yprob   = zeros(1, NoiseMcra.len);
% NoiseMcra.lamda_d   = zeros(1, NoiseMcra.len);
% NoiseMcra.b   = zeros(1, 2*NoiseMcra.w+1);    %hanning window coefficients

NoiseMcra.frm_cnt = NoiseMcra.frm_cnt + 1;

if (length(energy) ~= NoiseMcra.len)
    disp('error: length(energy) ~= NoiseMcra.len');
end
noise = zeros(1, NoiseMcra.len);
Sf = energy;

% Sf = zeros(1, NoiseMcra.len);
% for n = 1 : NoiseMcra.len
%     for j = 1 : 2*NoiseMcra.w+1
%         if ( n-NoiseMcra.w+j-1 > 0 && n-NoiseMcra.w+j-1 <= NoiseMcra.len )
%             Sf(n) = Sf(n) + NoiseMcra.b(j) * energy(n-NoiseMcra.w+j-1);
%         end
%     end
% end

if NoiseMcra.first == 1
    NoiseMcra.S = Sf;
    NoiseMcra.Smin = Sf;
    NoiseMcra.Stmp = Sf;
else
    NoiseMcra.S = 0.8 * NoiseMcra.S + 0.2 * Sf;
    NoiseMcra.Smin = min([NoiseMcra.Smin; NoiseMcra.S], [], 1);
    NoiseMcra.Stmp = min([NoiseMcra.Stmp; NoiseMcra.S], [], 1);
end

if NoiseMcra.frm_cnt == NoiseMcra.L
    NoiseMcra.Smin = min([NoiseMcra.Stmp; NoiseMcra.S], [], 1);
    NoiseMcra.Stmp = NoiseMcra.S;
    NoiseMcra.frm_cnt = 0;
end

Sr = NoiseMcra.S ./ NoiseMcra.Smin;
Yind = Sr > 2;

if NoiseMcra.first == 1
    NoiseMcra.Yprob = Yind;
else
    NoiseMcra.Yprob = 0.2 * NoiseMcra.Yprob + 0.8 * Yind;
end
alpha_td = NoiseMcra.alpha + (1-NoiseMcra.alpha) * NoiseMcra.Yprob;
if NoiseMcra.first == 1
    NoiseMcra.lamda_d = energy;
else
    NoiseMcra.lamda_d = alpha_td .* NoiseMcra.lamda_d + (1 - alpha_td) .* energy;
end
NoiseMcra.first = 0;
noise = NoiseMcra.lamda_d;