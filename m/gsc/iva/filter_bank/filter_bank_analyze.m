function X = filter_bank_analyze(x, frame_size, over_sample_ratio, win)

num_bands = frame_size * over_sample_ratio;
filter_size = length(win);
num_blocks = filter_size / num_bands;

if over_sample_ratio < 1, error('OvrSampRatio must be > 1'); end
if over_sample_ratio ~= 2, win = win/sqrt(over_sample_ratio/2); end

assert(2^round(log2(over_sample_ratio)) == over_sample_ratio);
assert(filter_size == length(win));

win = win(:);
% win = flipud(win);

% framing
x2 = buffer(x, filter_size, filter_size - frame_size);

% applying window funciton
x2 = x2 .* repmat(win, 1, size(x2, 2));

% folding
x3 = x2(1:num_bands, :);
for i=1:num_blocks-1
    x3 = x3 + x2(num_bands*i + (1:num_bands), :);
end

% fft
X = fft(x3);
% drop negative frequency component because of conjugate symmetry
X(num_bands/2+2:end, :) = []; 

end
