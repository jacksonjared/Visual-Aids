close all; clear; clc;

f_signal = 1e0;
f_carrier = 1e1;    % carrier frequency should be higher than signal frequency

Fs = 1e1 * f_carrier;   % Sampling frequency, must be at least double carrier frequency
L = 1e1 * Fs;                % number of sample periods, as this increases, fft becomes more accurate
t = [0:L-1] / Fs;       % Time vector

lower = 0;
upper = 2 / f_signal;

# information signal
signal = cos(2 * pi * f_signal * t);

# carrier signal
carrier = cos(2 * pi * f_carrier * t);

DSB = signal .* carrier;

AM = (1 + signal) .* carrier;

fDSB = fft(DSB);

fAM = fft(AM);

f = Fs * [0:L/2] / L;

P2 = abs(fDSB/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

subplot(2,2,1);
plot(t, carrier, 'g', t, DSB, 'b', t, signal, 'r');
xlim([lower upper]);
ylim([(min(DSB) - abs(min(DSB) * 0.1)) max(DSB) * 1.1])
title('Signal and DSB in time domain');
xlabel('time [s]');
ylabel('amplitude [V]');
legend('Carrier', 'DSB', 'Signal');

subplot(2,2,2);
plot(t, AM, 'g', t, signal, 'r');
xlim([lower upper]);
ylim([(min(AM) - abs(min(AM) * 0.1)) max(AM) * 1.1])
title('Signal and AM in time domain');
xlabel('time [s]');
ylabel('amplitude [V]');
legend('AM', 'Signal');

subplot(2,2,3);
plot(f, P1, 'r');
ylim([0 max(P1) * 1.1])
title('DSB in frequency domain');
xlabel('frequencu [hz]');
ylabel('one sided power spectrum');
legend('fDSB');

subplot(2,2,4);
plot(f, abs(fAM), 'r');
ylim([0 max(fAM) * 1.1])
title('AM in frequency domain');
xlabel('frequencu [hz]');
ylabel('amplitude');
legend('fAM');