% The basis Consists of polynomial functions shown as:
% p(x) = {x_1^a_1 × x_1^a_2 × ... × x_n^a_n | a_1 + ... a_n <= N}
clearvars, clc;
close all;
tf = 20;
nsteps = 1e3;
tspan = linspace(0, tf, nsteps);
initCond = rand*[pi/4 0 pi/3 0];
[T, X] = doublepend(tspan, initCond, false);

%% Initializing DMD [as per Chapter 16]
% first M measurements are used to boot the 'snapshot'
% psi = x itself in this example
M = 10;
X = X';
degree = 3;
X_dmd = X;
J = zeros(numel(tspan), 1);
for i = 1:numel(tspan)-M
    X_curr = X(:, i:i+M-1);
    X_next = X(:, i+1:i+M);
    [U, w_star] = koopID(X_curr, X_next, degree);
    
%     [z_tmp, ctrID] = monomialBasisEvaluator(X_curr(:,end), degree);
%     xx_tmp = U * z_tmp;
%     xx = xx_tmp(ctrID);
%     X_dmd = [X_dmd, xx];
%     evaltmp = monomialBasisEvaluator(X(:, i+M), degree);
%     J(i) = norm(xx_tmp - evaltmp);
    [z_tmp, ctrID] = monomialBasisEvaluator(X_curr(:,end), degree);
    xx_tmp = w_star * z_tmp;
    X_dmd(:,i+M) = xx_tmp;
     J(i) = norm(xx_tmp - X(:,i+M));
end

%% Plotting
figure, plot(T, X_dmd, 'LineWidth', 2.5);
hold on
plot(T, X,'--', 'LineWidth', 2.5);
grid minor;
grid on;
title('States in DMD v/s analytic [double pendulum]', 'FontSize', 16)
xlabel('Time [s]', 'FontSize', 14)
ylabel('Solutions [rad or rad/s]', 'FontSize', 14)
legend('\theta_2^{DMD}','d\theta_2^{DMD}/dt','\theta_1^{DMD}','d\theta_1^{DMD}/dt',...
    '\theta_2','d\theta_2/dt','\theta_1','d\theta_1/dt')