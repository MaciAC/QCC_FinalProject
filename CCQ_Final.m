syms theta

angle = theta/2;
Qs_0 = [1; 0];
Qs_1 = [0; 1];
theta_0 = [cos(angle); sin(angle)];
theta_1 = [cos(angle); -sin(angle)];

rho_XB = 1/2 * kron(Qs_0*Qs_0', theta_0*theta_0') + 1/2*kron(Qs_1*Qs_1', theta_1*theta_1');

rho_B = ptrace(rho_XB,1,[2,2]);

H_rho_B = Q_entropy_2D(rho_B);

I_XB = 1/2 * (H_rho_B - Q_entropy_2D(rho_B^0));

fplot(I_XB, [0, pi])

hold on

hold off
xticks([0, pi/2, pi])
xticklabels({'0', '\pi/2', '\pi'})
