%%% Just so we are all clear, the FFT is just a clever way to calculate the
%   DFT, it's not a seperate thing.

%%% Just MATLAB things %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear variables; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Sandbox %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Discrete, time domain signal
x = @(n) cos(2 .* pi .* 0.1 .* n) + cos(2 .* pi .* 0.2 .* n);

%%% Sampling frequency in Hz
Fs = 48e3;

%%% Sampling Period
Ts = 1 / Fs;

%%% Number of samples
N = 100;

%%% Padding length for FFT, only add a positive integer to N
N_padded = N + 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Time value of discrete samples
n = (0:N-1);

%%% Freq domain of the time signal
y = fft(x(n), N_padded);

figure;
subplot(2,1,1);
stem(n, x(n));
xlabel("n [sample]");
ylabel("Signal Amplitude");
title("Discrete Time Domain");

subplot(2,1,2);
%%% abs(y) to get the magnitude response. you have to multiply by 2 and
%   divide by N to normalize the magnitude response such that the value at
%   each frequency is the coefficient of that corresponding sine/cosine
%   wave. You have to multiply the x axis by 2*pi and divide by the length
%   to scale it to correspond to angular frequency.
stem((0:N_padded - 1) * 2 * pi / N_padded, abs(y) / N * 2, 'LineStyle', ...
                                                        '-', 'Color', 'r');
xlim([0 pi]);
xticks([0 pi/4 pi/2 (3*pi/4) pi]);
xticklabels(["0" "\pi/4" "\pi/2" "3\pi/4" "\pi"]);
xlabel("Normalized Angular Frequency [rad/s]");
ylabel("Frequency Component Amplitude");
title("Frequency Domain");

sgtitle("FS: " + num2str(Fs))