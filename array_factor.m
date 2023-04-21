close all; clear variables; clc;


N = 5;
psi = 0;
beta_d = 1.6*pi;


phi = linspace(0, 2*pi, 10000);

psi = psi * pi / 180;

chi = beta_d .* cos(phi) - psi;

i = 0:1:(N-1);
AF2 = (1/N) .* sum( exp(1j .* i' * chi) );

polarplot(phi, abs(AF2));