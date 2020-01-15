function [noise, NoiseMcra, Sr] = mcra(energy, NoiseMcra)

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

alpha_s = 0.8;
gamma = 0.998;
beta = 0.9;
alpha_p = 0.2;

mid = 96;
hi = 257;



if (length(energy) ~= NoiseMcra.len)
    disp('error: length(energy) ~= NoiseMcra.len');
end
Sf = energy;

if NoiseMcra.first == 1
    NoiseMcra.S = Sf;
    NoiseMcra.Smin = Sf;
else
    tmpS=NoiseMcra.S;
    % (9.61)
    NoiseMcra.S = alpha_s * NoiseMcra.S + (1 - alpha_s) * Sf;
    for n = 1 : NoiseMcra.len
        if(NoiseMcra.Smin(n) < NoiseMcra.S(n))
            % (9.25)
            NoiseMcra.Smin(n)=gamma*NoiseMcra.Smin(n)+(1 - gamma)/ (1 - beta) *(NoiseMcra.S(n)- beta *tmpS(n));
        else
            NoiseMcra.Smin(n) = NoiseMcra.S(n);
        end
    end
    
end



% (9.57)
Sr = NoiseMcra.S ./ NoiseMcra.Smin;

% (9.58) & %(9.60)
for n=1:mid
    if((Sr(n) > 2))
        NoiseMcra.P(n) = 1;
    else
        NoiseMcra.P(n) = 0;
    end
end
for n=mid + 1:hi
    if((Sr(n) > 5))
        NoiseMcra.P(n) = 1;
    else
        NoiseMcra.P(n) = 0;
    end
end

% (9.59)
NoiseMcra.Yprob = alpha_p * NoiseMcra.Yprob + (1 - alpha_p) *NoiseMcra.P;

% (9.54)
alpha_td = NoiseMcra.alpha + (1-NoiseMcra.alpha) * NoiseMcra.Yprob;

% (9.53)
NoiseMcra.lamda_d = alpha_td .* NoiseMcra.lamda_d + (1 - alpha_td) .* energy;

NoiseMcra.first = 0;
noise = NoiseMcra.lamda_d;