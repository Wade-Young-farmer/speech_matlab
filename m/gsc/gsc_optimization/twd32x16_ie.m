function [twd]=twd32x16_ie(N)
twd = exp(-2j*pi*[1;2;3]*(0:N/4-1)/N);
twd=twd.';
%size(twd)
twd = reshape([imag(twd(:).');real(twd(:).')],1,2*numel(twd)); 
twd = int16(round(pow2(twd,15)));

