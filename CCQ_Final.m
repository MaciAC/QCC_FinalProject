%% Plot I(X;B)
syms theta
angle = theta/2;
X_0 = [1; 0];
X_1 = [0; 1];
theta_0 = [cos(angle); sin(angle)];
theta_1 = [cos(angle); -sin(angle)];

rho_XB = 1/2 *( kron(X_0*X_0', theta_0*theta_0') + kron(X_1*X_1', theta_1*theta_1'));

rho_B = ptrace(rho_XB,1,[2,2]);

H_rho_B = Q_entropy_2D(rho_B);



I_XB = H_rho_B - 1/2*(Q_entropy_2D(theta_0*theta_0') + Q_entropy_2D(theta_1*theta_1'));

fplot(I_XB, [0, pi])

hold on

xticks([0, pi/2, pi])
xticklabels({'0', '\pi/2', '\pi'})

%% Plot I(X;Y)
P_e = 0.5*(1-sin(theta));
H_P = -P_e*log2(P_e) - (1-P_e)*log2(1-P_e);
I_Y = 1 - H_P;

fplot(I_Y, [0, pi])

hold on

hold off
legend('I(X;B)' ,'I(X;Y)')
xticks([0, pi/2, pi])
xticklabels({'0', '\pi/2', '\pi'})

%% POVM check

phi_0 = kron(kron(theta_0, theta_0),theta_0);
phi_1 = kron(kron(theta_0, theta_1),theta_1);
phi_2 = kron(kron(theta_1, theta_0),theta_1);
phi_3 = kron(kron(theta_1, theta_1),theta_1);
phis = [phi_0 phi_1 phi_2 phi_3];

P_x = 1/4;
vec_x = [[1; 0; 0; 0] [0; 1; 0; 0] [0; 0; 1; 0] [0; 0; 0; 1]]; 

rho_XB3 = zeros(32);
rho_B3 = zeros(8);
for i=[1:4]
    L_kron =  P_x*vec_x(:,i)*vec_x(:,i)';
    R_kron = phis(:,i)*phis(:,i)';
    rho_B3 = rho_B3 + R_kron;
    rho_XB3 = rho_XB3 + kron(L_kron, R_kron);
end

POVM_sum = zeros(8);
sqrt_rho_B3 = pinv(sqrt(rho_B3));
for i = [1:4]
    POVM_matrices(:,:,i) = 0.25 * sqrt_rho_B3*phis(:,i)*phis(:,i)';
    POVM_sum = POVM_sum + POVM_matrices(:,:,i);
end

%% Plot I(X;B3)
syms theta

angle = theta/2;
theta_0 = [cos(angle); sin(angle)];
theta_1 = [cos(angle); -sin(angle)];
phi_0 = kron(kron(theta_0, theta_0),theta_0);
phi_1 = kron(kron(theta_0, theta_1),theta_1);
phi_2 = kron(kron(theta_1, theta_0),theta_1);
phi_3 = kron(kron(theta_1, theta_1),theta_1);
phis = [phi_0 phi_1 phi_2 phi_3];

P_x = 1/4;
vec_x = [[1; 0; 0; 0] [0; 1; 0; 0] [0; 0; 1; 0] [0; 0; 0; 1]];

rho_XB3 = zeros(32);
for i = [1:4]
   rho_XB3 = rho_XB3 + kron(P_x*vec_x(:,i)*vec_x(:,i)', phis(:,i)*phis(:,i)');
end

rho_B3 = ptrace(rho_XB3, 1, [2,4,4]);

H_rho_B3 = Q_entropy_2D(rho_B3);

H_condition_XB3 = 0
for i =[1:4]
    H_condition_XB3 = H_condition_XB3 + P_x*Q_entropy_2D(phis(:,i)*phis(:,i)');
end
I_XB3 = H_rho_B3 - H_condition_XB3;

fplot(I_XB3, [0, pi])

hold on

xticks([0, pi/2, pi])
xticklabels({'0', '\pi/2', '\pi'})