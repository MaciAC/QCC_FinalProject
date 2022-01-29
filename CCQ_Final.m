%% Plot I(X;B)
syms theta

angle = theta/2;
Qs_0 = [1; 0];
Qs_1 = [0; 1];
theta_0 = [cos(angle); sin(angle)];
theta_1 = [cos(angle); -sin(angle)];

rho_XB = 1/2 * kron(Qs_0*Qs_0', theta_0*theta_0') + 1/2 * kron(Qs_1*Qs_1', theta_1*theta_1');

rho_B = ptrace(rho_XB,1,[2,2]);

H_rho_B = Q_entropy_2D(rho_B);



I_XB = H_rho_B - 1/2*(Q_entropy_2D(theta_0*theta_0') + Q_entropy_2D(theta_1*theta_1'));

fplot(I_XB, [0, pi])

hold on

xticks([0, pi/2, pi])
xticklabels({'0', '\pi/2', '\pi'})

%% Plot I(X;Y)
P_e = 0.5*(1-sin(theta));
H_P = -P_e*log(P_e) - (1-P_e)*log(1-P_e);
I_Y = 1 - H_P;

fplot(I_Y, [0, pi])

hold on

hold off
xticks([0, pi/2, pi])
xticklabels({'0', '\pi/2', '\pi'})