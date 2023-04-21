%%% Just so we are all clear, the FFT is just a clever way to calculate the
%   DFT, it's not a seperate thing.

%%% Just MATLAB things %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear variables; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Sandbox %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Time domain signal
x = @(t) cos(2 * pi * 4800 * t) + cos(2 * pi * 9600 * t);

%%% Sampling time in seconds
T = 0.1;

%%% Sampling frequency in Hz
Fs = 50e3;

%%% Sampling Period
Ts = 1 / Fs;

%%% Number of samples. This can either be set explicitly, or calculated
%   based on the sampling freq, Fs, and the sampling time, T
N = Fs * T;
% N = 1e6;

%%% Padding length for FFT
N_padded = N + 0;

%%% Set to 1 to normalize the frequency axis.
freq_norm = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Time value of discrete samples
time = (0:N-1) * Ts;

%%% Freq domain of the time signal
y = fft(x(time), N_padded);

%%% Freq axis scaling factor
freq_scale = (1 - freq_norm) * (Fs / N_padded) + freq_norm * ...
                                                            (1 / N_padded);
freq_axis = (1 - freq_norm) * (Fs / 2) + freq_norm * 0.5;

figure;
subplot(2,1,1);
plot(time, x(time));
xlabel("Time [s]");
ylabel("Signal Amplitude");
title("Time Domain");

subplot(2,1,2);
%%% abs(y) to get the magnitude response. you have to multiply by 2 and
%   divide by N to normalize the magnitude response such that the value at
%   each frequency is the coefficient of that corresponding sine/cosine
%   wave. You have to multiply the x axis by Fs and divide my the length to
%   scale it to correspond to frequency.
plot((0:N_padded - 1) * freq_scale, abs(y) / N * 2, 'LineStyle', '-', ...
                                                             'Color', 'r');
xlim([0 freq_axis]);
if (freq_norm)
    xlabel("Normalized Frequency [Hz]");
else
    xlabel("Frequency [Hz]");
end
ylabel("Frequency Component Amplitude");
title("Frequency Domain");

sgtitle("FS: " + num2str(Fs))