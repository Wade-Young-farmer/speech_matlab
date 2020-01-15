function [noise, NoiseMcra] = mcra(energy, NoiseMcra)

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

if (length(energy) ~= NoiseMcra.len)
    disp('error: length(energy) ~= NoiseMcra.len');
end
Sf = energy;

tmpS=NoiseMcra.S;
NoiseMcra.S = 0.8 * NoiseMcra.S + 0.2 * Sf;
    for n = 1 : NoiseMcra.len
        if(NoiseMcra.Smin(n) < NoiseMcra.S(n))
            NoiseMcra.Smin(n)=0.998*NoiseMcra.Smin(n)+0.02*(NoiseMcra.S(n)-0.9*tmpS(n));
        else
            NoiseMcra.Smin(n) = NoiseMcra.S(n);
        end
    end


Sr = NoiseMcra.S ./ NoiseMcra.Smin;
for n=1:48
    if((Sr(n) > 2))
        NoiseMcra.P(n) = 1;
    else
        NoiseMcra.P(n) = 0;
    end
end
for n=49:129
    if((Sr(n) > 5))
        NoiseMcra.P(n) = 1;
    else
        NoiseMcra.P(n) = 0;
    end
end

NoiseMcra.Yprob = 0.2 * NoiseMcra.Yprob + 0.8 *NoiseMcra.P;

alpha_td = NoiseMcra.alpha + (1-NoiseMcra.alpha) * NoiseMcra.Yprob;
NoiseMcra.lamda_d = alpha_td .* NoiseMcra.lamda_d + (1 - alpha_td) .* energy;

NoiseMcra.first = 0;
noise = NoiseMcra.lamda_d;