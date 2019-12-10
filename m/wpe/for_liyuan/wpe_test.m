close all; clear; clc;
dbstop if error;

addpath('./filter_bank');

file_path = './wav_sample';
in_file = fullfile(file_path, 'sample_4ch.wav');
reverb_time = 500;

out_filename = ['wpe_out_rt' int2str(reverb_time) '.wav'];
out_file = fullfile(file_path, out_filename);

[x, fs] = audioread(in_file);
mic_num = size(x, 2);

frame_len = 10;
frame_size = fs * frame_len / 1000;
over_sample_ratio = 2;
subband_num = frame_size * over_sample_ratio / 2 + 1;
frame_num = floor(size(x, 1) / frame_size);

delay = 1;
h_size = frame_size * over_sample_ratio * 3;
h = filter_bank_win(h_size, delay);

X = zeros(mic_num, subband_num, frame_num);
for ii = 1:mic_num
    X(ii, :, :) = filter_bank_analyze(x(:, ii), frame_size, over_sample_ratio, h);
end

tic
Y = wpe_process(X, frame_len, reverb_time);
toc

Y = squeeze(Y(1, :, :));
y = filter_bank_synthesis(Y, frame_size, over_sample_ratio, h);

audiowrite(out_file, y, fs);

figure(1);clf
plot((0:length(x)-1)/fs, x(:, 1), 'b');
hold on
plot((0:length(y)-1)/fs, y, 'r');
grid on
