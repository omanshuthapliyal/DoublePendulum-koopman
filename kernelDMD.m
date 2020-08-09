function [Phi_x1 ,omega ,lambda ,Xi,Xkerneldmd] = kernelDMD(X1,X2,r,dt)
[n,m]=size(X1);
lambda = 0.5;
for j=1:m
    for jj=1:m
        A(j,jj)=exp(-lambda*abs((X1(:,jj)-X2(:,j)).'*(X1(:,jj)-X2(:,j)))^2);
        G(j,jj)=exp(-lambda*abs((X1(:,jj)-X1(:,j)).'*(X1(:,jj)-X1(:,j)))^2);
    end
end
    
[Q,Sig,~] = svd(G);
r = min(r, size(Q,2));
% Q_r = Q(:, 1:r); % truncate to rank-r
% S_r = Sig(1:r, 1:r);
% % V_r = V(:, 1:r);
% S_r = S_r.^(0.5);
% tempMat = pinv(S_r)*Q_r';
Sig = Sig.^(0.5);
tempMat = Sig \ Q';
K_hat = tempMat * A * tempMat';
[V,D] = eig(K_hat);
lambda = diag(D); % discrete -time eigenvalues
omega = log(lambda)/dt; % continuous-time eigenvalues
Phi_x1 = G * tempMat' * V;  % eigenfunctions evaluated at X1
Phi_x2 = A * tempMat' * V;  % eigenfunctions evaluated at X2

Xi = Phi_x1\X1';    % Xi has shape [xi_1,xi_2,...,xi_M]^T
% Xi = inv(V)*tempMat*X1';
Xi = Xi';   
%% DMD prediction for next time
Xkerneldmd = zeros(n,1);
time_dynamics = zeros(n, m);
k = size(Xi,2);
phi_x = Phi_x1(end,:);
t = (0:m -1)*dt; % time vector
for i = 1:k
    Xkerneldmd = Xkerneldmd + (Xi(:,i)./norm(Xi(:,i)) ).*exp(omega(i)*dt).*phi_x(i);
end
% for iter = 1:m
%     time_dynamics (:,iter )=(Xi.*exp(omega*t(iter )));
% end
% % Xdmd = Phi_x1 * time_dynamics ;
% Xdmd = time_dynamics * Phi_x1;
end