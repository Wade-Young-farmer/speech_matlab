function x = filter_bank_synthesis(X, frame_size, over_sample_ratio, win)

num_bands = frame_size * over_sample_ratio;
filter_size = length(win);
num_blocks = filter_size / num_bands;

if over_sample_ratio < 1, error('OvrSampRatio must be > 1'); end
if over_sample_ratio ~= 2, win = win/sqrt(over_sample_ratio/2); end

assert(2^round(log2(over_sample_ratio)) == over_sample_ratio);
assert(filter_size == length(win));

win = win(:);

% recover negative frequencies
X = [X; conj(flipud(X(2:num_bands/2, :)))];
[~, num_frames] = size(X);

% ifft
x2 = real(ifft(X));

% scale up because the scaling is included in the window function
x2 = x2 * num_bands;

% prepare output buffer
x = zeros(frame_size*(num_frames+num_blocks*over_sample_ratio-1),1);

% unfolding
x2 = repmat(x2, num_blocks, 1);

% apply window
x2 = x2.*repmat(win, 1, num_frames);

% overlap and add
for i=1:num_frames
    start_index = (i-1)*frame_size;
    x(start_index + (1:filter_size)) = x(start_index + (1:filter_size)) + x2(:,i);
end

end
