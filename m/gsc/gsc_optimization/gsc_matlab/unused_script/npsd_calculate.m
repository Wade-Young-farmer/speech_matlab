function [npsd, spsd, alpha] = npsd_calculate(input_spsd, ctrl_snr, m_ALPHAc)
len = 65;
ctrl_snr = zeros(len, 1);
spsd = zeros(len, 1);
input_spsd = zeros(len, 1);

% input_spsd
%% Calculate alpha_c
ALPHAcs = 0;
temp = 0;
for i = 1:len
    temp = temp + ctrl_snr(i);
    ALPHAcs = ALPHAcs + input_spsd(i);
end

ALPHAcs = ALPHAcs / (temp + 1.0000e-09);
ALPHAcs = ALPHAcs - 1;
ALPHAcs = ALPHAcs * ALPHAcs;
ALPHAcs = 1 + ALPHAcs;
ALPHAcs = 1 / ALPHAcs;

m_ALPHAc = 0.7 * m_ALPHAc + 0.3 * max(ALPHAcs, 0.7);

%% Calculate alpha_p
npsd = 0;
alpha = 0;
