function [twd]=twdf_ie(N)
twd = exp(-2j*pi*[1;2;3]*(0:N/4-1)/N);
size(twd)
twd = reshape([real(twd(:).');imag(twd(:).')],1,2*numel(twd));