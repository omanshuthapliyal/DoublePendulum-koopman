% The basis Consists of polynomial functions shown as:
% p(x) = {x_1^a_1 × x_1^a_2 × ... × x_n^a_n | a_1 + ... a_n <= N}
clearvars, clc;
close all;
tf = 20;
nsteps = 1e3;
tspan = linspace(0, tf, nsteps);
dt = tspan(2)-tspan(1);
initCond = rand*[pi/4 0 pi/3 0];
[T, X] = doublepend(tspan, initCond, false);
n = size(X,2);
M = 2;
X = X';
Psi_X = zeros();
lambda = 1;
X_dmd = X;
X_kerneldmd = X;
for i = 1:numel(tspan)-M
    X1 = X(:, i:i+M-1);
    X2 = X(:, i+1:i+M);
    [Phi ,omega ,lambda ,b,Xdmd] = DMD(X1,X2,10,dt);
    xx_dmd = Xdmd(:,end);
    X_dmd(:,i+M) = xx_dmd;
    
    [Phi ,omega ,lambda ,b,Xkerneldmd] = kernelDMD(X1,X2,4,dt);
    xx_kerneldmd = Xkerneldmd(:,end);
    X_kerneldmd(:,i+M) = xx_kerneldmd;
end

%% Plotting
figure, plot(T, X_dmd, ':', 'LineWidth', 2.5);
hold on
plot(T, X_kerneldmd,'-.', 'LineWidth', 2.5);
% plot(T, X,'--', 'LineWidth', 2.5);
grid minor;
grid on;
title('States in DMD v/s analytic [double pendulum]', 'FontSize', 16)
xlabel('Time [s]', 'FontSize', 14)
ylabel('Solutions [rad or rad/s]', 'FontSize', 14)
legend('\theta_2^{DMD}','d\theta_2^{DMD}/dt','\theta_1^{DMD}','d\theta_1^{DMD}/dt',...
    '\theta_2^{kernel}','d\theta_2^{kernel}/dt','\theta_1^{kernel}','d\theta_1^{kernel}/dt'); %,...
%     '\theta_2','d\theta_2/dt','\theta_1','d\theta_1/dt')