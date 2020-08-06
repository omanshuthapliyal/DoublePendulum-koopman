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
M = 10;
X = X';
Psi_X = zeros();
lambda = 1;
X_kerneldmd = X;
%% Gaussina kernel based EDMD 
for i = 1:numel(tspan)-M
    x = X(:, i:i+M-1);
    y = X(:, i+1:i+M);
    % A = G^dagger H
    % G = \psi_X^T \psi_X
    % H = \psi_X^T \psi_Y
    A = zeros(M,M);
    G = A;
    for ii = 1:M
        for jj = 1:M
            A(ii,jj) = gaussianKernel(x(:,jj),y(:,ii), lambda);
            G(ii,jj) = gaussianKernel(x(:,jj),x(:,ii), lambda);
        end
    end
    [V,Sig] = eig(G);
    K = (Sig*V') * A * (V/Sig);
    [U,S,V] = svd(x);
    r = size(U,2);
    U_r = U(:, 1:r); % truncate to rank-r
    S_r = S(1:r, 1:r);
    V_r = V(:, 1:r);
    [W,D] = eig(K);
    
%     Phi = y * V_r / S_r * W;
%     lambda = diag(D);
%     omega = log(lambda)/dt;
    
    x_tmp = x(:,end);

%     
    [Q,S,Z] = svd(G);
    S = S.^(0.5);
    K = (pinv(S)*Q') * A * (Q*pinv(S));
    x_tmp = x(:,end);
    [W,D,V] = eig(K);
    for i_p = 1:M
        tmp = [];
        for j_p=1:M
        	tmp = [tmp, gaussianKernel(x_tmp,x(:,j_p), lambda)];
        end
    end
    phi = tmp*Q*pinv(S)*V;
    
%     L = logm(K)/dt;
%     xx = L*X(:, i:i+M-1)';
    b = phi\x_tmp;
    xx = Phi*b.*exp(omega*tspan(i))
    xx = xx(end,:)';
    X_kerneldmd(:,i+M) = xx;
end
%% Initializing DMD [as per Chapter 16]
% first M measurements are used to boot the 'snapshot'
% psi = x itself in this example
% X = X';
degree = 3;
X_dmd = X;
J = zeros(numel(tspan), 1);
for i = 1:numel(tspan)-M
    X_curr = X(:, i:i+M-1);
    X_next = X(:, i+1:i+M);
    [Q, w_star] = koopID(X_curr, X_next, degree);
    
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
% figure, plot(T, X_dmd, ':', 'LineWidth', 2.5);
% hold on
plot(T, X_kerneldmd, '-', 'LineWidth', 2.5);
hold on
plot(T, X, ':','LineWidth', 2.5);
grid minor;
grid on;
title('States in DMD v/s analytic [double pendulum]', 'FontSize', 16)
xlabel('Time [s]', 'FontSize', 14)
ylabel('Solutions [rad or rad/s]', 'FontSize', 14)
legend('\theta_2^{DMD}','d\theta_2^{DMD}/dt','\theta_1^{DMD}','d\theta_1^{DMD}/dt',...
    '\theta_2','d\theta_2/dt','\theta_1','d\theta_1/dt')