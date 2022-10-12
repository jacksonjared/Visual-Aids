close all; clear variables; clc;

%%% vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
lower = 2;

upper = 4;
%%% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


f = logspace(lower, upper, 10000);
w = 2 .* pi .* f;
Wc = 10;
Sa = 1i .* w;
Sb = Sa ./ Wc;

%%% vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
voltageTransfer = ( 1i .* w .* 0.328 ./ 5.2e-3 ) ./ ( -(w).^2 + 1i .* w .* 0.328 ./ 5.2e-3 + (1 / (5.2e-3 * 10e-6)) );
%%% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


w = 2 .* pi .* f;
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