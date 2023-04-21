close all; clear variables; clc;

%%% vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
lower = -10;

upper = 0;
%%% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

f = logspace(lower, upper, 10000);
w = 2 .* pi .* f;

%%% vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
voltageTransfer = ( -10 ./ (1j .* w + 0.01) ) + ( (-10 .* (1j .* w).^2 .* w) ./ ( (1j .* w) .* ((1j .* w).^2 + (w).^2) .* (1j .* w + 0.01) ) )
%%% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

bodeplot = 20*log10(abs(voltageTransfer));

% subplot(2,2,1);
semilogx(f, bodeplot, 'r', f, -3.0103 + 0.*f, '--b');
xlim([10^(lower) 10^(upper)]);
ylim([min(bodeplot) (max(bodeplot) + 2)]);
grid on;
xlabel("Frequency [Hz]");
ylabel("Magnitude [dB]");

% subplot(2,2,2);
% semilogx(w, bodeplot, 'r', w, -3.0103 + 0.*w, '--b');
% xlim([10^(lower) 10^(upper)]);
% ylim([-200 1]);
% grid on;
% xlabel("Angular Frequency [rad/s]");
% ylabel("Magnitude [dB]");
% 
% subplot(2,2,[3,4]);
% semilogx(f, angle(voltageTransfer) .* 180 ./ pi, 'r');
% xlim([10^(lower) 10^(upper)]);
% ylim([-180 180]);
% grid on;
% xlabel("Frequency [dB]");
% ylabel("Phase angle in degrees");