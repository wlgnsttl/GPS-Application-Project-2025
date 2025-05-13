function[SatPos, dRel] = getSatPos_lab(eph, ieph, gs)
% eph, ieph, gs 받아 위성pos 반환 [xyz]

t_oe = eph(ieph, 01);
sqrtA = eph(ieph, 10);
e = eph(ieph, 11);
i_0 = eph(ieph, 12);
omega = eph(ieph, 13);
Omega_0 = eph(ieph, 14);
M_0 = eph(ieph, 15);
i_dot = eph(ieph, 16);
Omega_dot = eph(ieph, 17);
dn = eph(ieph, 18);
C_uc = eph(ieph, 20);
C_us = eph(ieph, 21);
C_rc = eph(ieph, 22);
C_rs = eph(ieph, 23);
C_ic = eph(ieph, 24);
C_is = eph(ieph, 25);

mu = 3.986005e14;
Omega_dot_e = 7.2921151467e-5;

a = sqrtA^2;
n_0 = sqrt(mu/a^3);
n = n_0 + dn;
t_k = gs - t_oe;
M_k = M_0 + n*t_k;
E_k = solveKeplerEq(M_k, e);
f_k =  atan2((sqrt(1 - e^2)*sin(E_k)), (cos(E_k) - e));

phi_k = f_k + omega;

delta_uk = C_us*sin(2*phi_k) + C_uc*cos(2*phi_k);
delta_rk = C_rs*sin(2*phi_k) + C_rc*cos(2*phi_k);
delta_ik = C_is*sin(2*phi_k) + C_ic*cos(2*phi_k);

u_k = phi_k + delta_uk;
r_k = a*(1 - e*cos(E_k)) + delta_rk;
i_k = i_0 + delta_ik + i_dot*t_k;

x_kp = r_k*cos(u_k);
y_kp = r_k*sin(u_k);

Omega_k = Omega_0 + (Omega_dot - Omega_dot_e)*t_k - Omega_dot_e*t_oe;

x_k = x_kp*cos(Omega_k) - y_kp*cos(i_k)*sin(Omega_k);
y_k = x_kp*sin(Omega_k) + y_kp*cos(i_k)*cos(Omega_k);
z_k = y_kp*sin(i_k);


SatPos = [x_k, y_k, z_k];

F = -4.442807633e-10;
dRel = F*e*sqrtA*sin(E_k);
