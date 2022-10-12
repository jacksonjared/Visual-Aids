%%% This program attempts to plot the vector field associated with a transverse
%%% magnetic wave travelling down a rectangular waveguide in free space at time
%%% 0.

%%% Note: E field and H field magnitudes are correct relative to themselves, but
%%% not to each other.

%%% TODO:
%%% Make X, Y, Z resolution independent

%% Just MATLAB things %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear variables; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Sandbox %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% The mode of the wave in the x-axis
m = 1;

%%% The mode of the wave in the y-axis
n = 1;

%%% Length of waveguide in x-axis
a = 1;

%%% Length of waveguide in y-axis
b = 1;

%%% The plane to plot in wavelengths [0 1], if this number is negative, it
%%% will step through all planes
plane = -0.25;

%%% Angular freq of the wave
omega = 4e10;

%%% Spacial resolution
resS = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Vacuum permittivity
mu0 = pi * 4e-7;

%%% Vacuum permeability
epsilon0 = 8.8541878128e-12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Maths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Helper constants
M = m * pi / a;
N = n * pi / b;

%%% Beta of the wave if it was a plane wave
beta0 = omega * sqrt(mu0 * epsilon0);

%%% The beta of the wave as it is in the waveguide
betaMN = sqrt( (omega^2 * mu0 * epsilon0) - M^2 - N^2 );

%%% The wavelength of the wave as it is in the waveguide
lambaMN = 2 * pi / betaMN;

%%% The lowest angular frequency that can travel in the given waveguide
omega_min = (1 / sqrt(mu0 * epsilon0)) * sqrt(M^2 + N^2);

%%% Set up the grid of points to plot at based on the waveguide and spacial
%%% resolution
[Z, X, Y] = meshgrid(linspace(0,lambaMN, resS+1), linspace(0,a,resS+1), ...
                                                          linspace(0,b,resS+1));

%%% What index that plane is closest to
index = round(plane * resS) + 1;

%%% Pre define the z axis term since it is the same for every component
Kt = exp(-1j .* betaMN .* Z);

%%% Define a constant to save calculation time later
Ke = (-1j * betaMN * pi) / ( (omega^2 * mu0 * epsilon0) - betaMN^2 );
Kh = (1j * omega * epsilon0 * pi) /( (omega^2 * mu0 * epsilon0) - betaMN^2 );

%%% Z-axis component of H field (0 because it is a TM wave)
Hz = zeros(resS+1,resS+1,resS+1);

%%% Y-axis component of H field
Hy = (-Kh * m / a) .* cos(M .* X) .* sin(N .* Y) .* Kt;

%%% X-axis component of H field
Hx = (Kh * n / b) .* sin(M .* X) .* cos(N .* Y) .* Kt;

%%% Z axis component of E field
Ez = sin(M .* X) .* sin(N .* Y) .* Kt;

%%% Y axis component of E field
Ey = (Ke * n / b) .* sin(M .* X) .* cos(N .* Y) .* Kt;

%%% X axis component of E field
Ex = (Ke * m / a) .* cos(M .* X) .* sin(N .* Y) .* Kt;

%%% Scale H to avoid overlapping, while keeping the correct relative amplitude
H_scale = min([a/resS b/resS]) / sqrt(2);
H_max = max( [max(real(Hz),[],'all') max(real(Hx),[],'all') max(real(Hy),[], ...
                                                            'all')] ) / H_scale;

%%% Scale E to avoid overlapping, while keeping the correct relative amplitude
E_scale = min([a/resS b/resS]) / sqrt(3);
E_max = max( [max(real(Ez),[],'all') max(real(Ex),[],'all') max(real(Ey),[], ...
                                                            'all')] ) / E_scale;

%%% Just plot the scaled, real part of H
plot_2d_Hy = real(Hy) ./ H_max;
plot_2d_Hx = real(Hx) ./ H_max;

%%% Just plot the scaled, real part of E
plot_2d_Ez = real(Ez) ./ E_max;
plot_2d_Ey = real(Ey) ./ E_max;
plot_2d_Ex = real(Ex) ./ E_max;

%%% Increase limits to see vectors
a10 = 0.1 * a;
b10 = 0.1 * b;
L10 = 0.1 * lambaMN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure();

%%% XY plane Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,1);

%%% Allows for greek characters in labels
set(gca, 'TickLabelInterpreter', 'latex');

%%% Configure x-axis
xlabel('Z-axis');
xticks([0 lambaMN/4 lambaMN/2 3*lambaMN/4 lambaMN])
xticklabels({'$0$','$\frac{\lambda}{4}$','$\frac{\lambda}{2}$', ...
                                            '$\frac{3\lambda}{4}$','$\lambda$'})

%%% Configure y-axis
ylabel('X-axis');
yticks([0 a/4 a/2 3*a/4 a])
yticklabels({'$0$','$\frac{a}{4}$','$\frac{a}{2}$','$\frac{3a}{4}$','$a$'})

%%% Configure z-axis
zlabel('Y-axis');
zticks([0 b/4 b/2 3*b/4 b])
zticklabels({'$0$','$\frac{b}{4}$','$\frac{b}{2}$','$\frac{3b}{4}$','$b$'})

%%% Set limits and viewport
axis square;
grid on;
axis([(-L10) (lambaMN+L10) (-a10) (a+a10) (-b10) (b+b10)]);
view([-1, -1, 1]);

hold on;

subplot(2,3,4)

%%% Allows for greek characters in labels
set(gca, 'TickLabelInterpreter', 'latex');

%%% Configure x-axis
xlabel('X-axis');
xticks([0 a/4 a/2 3*a/4 a])
xticklabels({'$0$','$\frac{a}{4}$','$\frac{a}{2}$','$\frac{3a}{4}$','$a$'})
set(gca, "XDir", 'reverse');

%%% Configure y-axis
ylabel('Y-axis');
yticks([0 b/4 b/2 3*b/4 b])
yticklabels({'$0$','$\frac{b}{4}$','$\frac{b}{2}$','$\frac{3b}{4}$','$b$'})

%%% Set limits and viewport
axis square;
grid on;
axis([(-a10) (a+a10) (-b10) (b+b10)]);
view(2);

hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% YZ plane Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,2);

%%% Allows for greek characters in labels
set(gca, 'TickLabelInterpreter', 'latex');

%%% Configure x-axis
xlabel('Z-axis');
xticks([0 lambaMN/4 lambaMN/2 3*lambaMN/4 lambaMN])
xticklabels({'$0$','$\frac{\lambda}{4}$','$\frac{\lambda}{2}$', ...
                                            '$\frac{3\lambda}{4}$','$\lambda$'})

%%% Configure y-axis
ylabel('X-axis');
yticks([0 a/4 a/2 3*a/4 a])
yticklabels({'$0$','$\frac{a}{4}$','$\frac{a}{2}$','$\frac{3a}{4}$','$a$'})

%%% Configure z-axis
zlabel('Y-axis');
zticks([0 b/4 b/2 3*b/4 b])
zticklabels({'$0$','$\frac{b}{4}$','$\frac{b}{2}$','$\frac{3b}{4}$','$b$'})

%%% Set limits and viewport
axis square;
grid on;
axis([(-L10) (lambaMN+L10) (-a10) (a+a10) (-b10) (b+b10)]);
view([-1, -1, 1]);

hold on;

subplot(2,3,5)

%%% Allows for greek characters in labels
set(gca, 'TickLabelInterpreter', 'latex');

%%% Configure x-axis
xlabel('Z-axis');
xticks([0 lambaMN/4 lambaMN/2 3*lambaMN/4 lambaMN])
xticklabels({'$0$','$\frac{\lambda}{4}$','$\frac{\lambda}{2}$', ...
                                            '$\frac{3\lambda}{4}$','$\lambda$'})

%%% Configure y-axis
ylabel('Y-axis');
yticks([0 b/4 b/2 3*b/4 b])
yticklabels({'$0$','$\frac{b}{4}$','$\frac{b}{2}$','$\frac{3b}{4}$','$b$'})

%%% Set limits and viewport
axis square;
grid on;
axis([(-L10) (lambaMN+L10) (-b10) (b+b10)]);
view(2);

hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% XZ plane Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,3);

%%% Allows for greek characters in labels
set(gca, 'TickLabelInterpreter', 'latex');

%%% Configure x-axis
xlabel('Z-axis');
xticks([0 lambaMN/4 lambaMN/2 3*lambaMN/4 lambaMN])
xticklabels({'$0$','$\frac{\lambda}{4}$','$\frac{\lambda}{2}$', ...
                                            '$\frac{3\lambda}{4}$','$\lambda$'})

%%% Configure y-axis
ylabel('X-axis');
yticks([0 a/4 a/2 3*a/4 a])
yticklabels({'$0$','$\frac{a}{4}$','$\frac{a}{2}$','$\frac{3a}{4}$','$a$'})

%%% Configure z-axis
zlabel('Y-axis');
zticks([0 b/4 b/2 3*b/4 b])
zticklabels({'$0$','$\frac{b}{4}$','$\frac{b}{2}$','$\frac{3b}{4}$','$b$'})

%%% Set limits and viewport
axis square;
grid on;
axis([(-L10) (lambaMN+L10) (-a10) (a+a10) (-b10) (b+b10)]);
view([-1, -1, 1]);

hold on;

subplot(2,3,6)

%%% Allows for greek characters in labels
set(gca, 'TickLabelInterpreter', 'latex');

%%% Configure x-axis
xlabel('Z-axis');
xticks([0 lambaMN/4 lambaMN/2 3*lambaMN/4 lambaMN])
xticklabels({'$0$','$\frac{\lambda}{4}$','$\frac{\lambda}{2}$', ...
                                            '$\frac{3\lambda}{4}$','$\lambda$'})

%%% Configure y-axis
ylabel('X-axis');
yticks([0 a/4 a/2 3*a/4 a])
yticklabels({'$0$','$\frac{a}{4}$','$\frac{a}{2}$','$\frac{3a}{4}$','$a$'})

%%% Set limits and viewport
axis square;
grid on;
axis([(-L10) (lambaMN+L10) (-a10) (a+a10)]);
view(2);

hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (plane < 0)
  start = 0;
  stop = inf;
else
  start = index;
  stop = index;
end

which_3d = 2;

for loop = start:stop

  index = mod(loop, resS + 1) + 1;

  if(index == 1)
      which_3d = mod(which_3d + 1, 3);
  end

%%% 3D Overview %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subplot(2,3,1)
    %%% E field XY plotting
    E_arrows_1 = quiver3(Z(:,index,:), X(:,index,:), Y(:,index,:), ...
                         plot_2d_Ez(:,index,:), plot_2d_Ex(:,index,:), ...
                         plot_2d_Ey(:,index,:),'b', 'AutoScale', 'off', ...
                                                        'ShowArrowHead', 'off');

    %%% H field XY plotting
    H_arrows_1 = quiver3(Z(:,index,:), X(:,index,:), Y(:,index,:), ...
                         Hz(:,index,:), plot_2d_Hx(:,index,:), ...
                         plot_2d_Hy(:,index,:), 'r', 'AutoScale', 'off', ...
                                                        'ShowArrowHead', 'off');

  subplot(2,3,2)
    %%% E field YZ plotting
    E_arrows_2 = quiver3(Z(index,:,:), X(index,:,:), Y(index,:,:), ...
                         plot_2d_Ez(index,:,:), plot_2d_Ex(index,:,:), ...
                         plot_2d_Ey(index,:,:),'b', 'AutoScale', 'off', ...
                                                        'ShowArrowHead', 'off');

    %%% H field YZ plotting
    H_arrows_2 = quiver3(Z(index,:,:), X(index,:,:), Y(index,:,:), ...
                         Hz(index,:,:), plot_2d_Hx(index,:,:), ...
                         plot_2d_Hy(index,:,:), 'r', 'AutoScale', 'off', ...
                                                        'ShowArrowHead', 'off');

  subplot(2,3,3)
    %%% E field XZ plotting
    E_arrows_3 = quiver3(Z(:,:,index), X(:,:,index), Y(:,:,index), ...
                         plot_2d_Ez(:,:,index), plot_2d_Ex(:,:,index), ...
                         plot_2d_Ey(:,:,index),'b', 'AutoScale', 'off', ...
                                                        'ShowArrowHead', 'off');

    %%% H field XZ plotting
    H_arrows_3 = quiver3(Z(:,:,index), X(:,:,index), Y(:,:,index), ...
                         Hz(:,:,index), plot_2d_Hx(:,:,index), ...
                         plot_2d_Hy(:,:,index), 'r', 'AutoScale', 'off', ...
                                                        'ShowArrowHead', 'off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% XY Plane %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subplot(2,3,4)

  %%% E field plotting
  E_arrows_XY = quiver(X(:,index,:), Y(:,index,:), plot_2d_Ex(:,index,:), ...
         plot_2d_Ey(:,index,:),'b', 'AutoScale', 'off', 'ShowArrowHead', ...
                                                                          'on');

  %%% H field plotting
  H_arrows_XY = quiver(X(:,index,:), Y(:,index,:), plot_2d_Hx(:,index,:), ...
         plot_2d_Hy(:,index,:),'r', 'AutoScale', 'off', 'ShowArrowHead', ...
                                                                          'on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% YZ Plane %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subplot(2,3,5)

  %%% E field plotting
  E_arrows_YZ = quiver(Z(index,:,:), Y(index,:,:), plot_2d_Ez(index,:,:), ...
         plot_2d_Ey(index,:,:),'b', 'AutoScale', 'off', 'ShowArrowHead', ...
                                                     'on', 'MaxHeadSize', 0.01);

  %%% H field plotting
  H_arrows_YZ = quiver(Z(index,:,:), Y(index,:,:), Hz(index,:,:), ...
         plot_2d_Hy(index,:,:),'r', 'AutoScale', 'off', 'ShowArrowHead', ...
                                                     'on', 'MaxHeadSize', 0.01);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% XZ Plane %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subplot(2,3,6)

  %%% E field plotting
  E_arrows_XZ = quiver(Z(:,:,index), X(:,:,index), plot_2d_Ez(:,:,index), ...
         plot_2d_Ex(:,:,index),'b', 'AutoScale', 'off', 'ShowArrowHead', ...
                                                     'on', 'MaxHeadSize', 0.01);

  %%% H field plotting
  H_arrows_XZ = quiver(Z(:,:,index), X(:,:,index), Hz(:,:,index), ...
         plot_2d_Hx(:,:,index),'r', 'AutoScale', 'off', 'ShowArrowHead', ...
                                                     'on', 'MaxHeadSize', 0.01);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%% Prevents looping if on a single plane is desired
  if(plane > 0)
    break;
  end

  pause(0.1);

  %%% Clean up the graphs
  delete(E_arrows_1);
  delete(H_arrows_1);
  delete(E_arrows_2);
  delete(H_arrows_2);
  delete(E_arrows_3);
  delete(H_arrows_3);
  delete(E_arrows_XY);
  delete(H_arrows_XY);
  delete(E_arrows_YZ);
  delete(H_arrows_YZ);
  delete(E_arrows_XZ);
  delete(H_arrows_XZ);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%