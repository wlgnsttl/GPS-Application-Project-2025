function SatVel = getSatVel(eph, ieph, gs)
%--------------------------------------------------------------
%  eph   : RINEX nav 행렬 (각·각속도 = rad / rad/s)
%  ieph  : 사용할 행(index)
%  gs    : 계산 epoch (GPS Time, s)
%  출력  : SatVel = [vx vy vz]  (m/s)
%  단위  : 모든 각 = rad, 각속도 = rad/s
%--------------------------------------------------------------

% --- 상수
mu = 3.986005e14;            % 지구 중력상수 (m^3/s^2)
we = 7.2921151467e-5;        % 지구 자전속도 (rad/s)

% --- 1. 방송 파라미터(라디안 그대로) --------------------------
t_oe   = eph(ieph,  1);       % reference epoch
sqrt_a = eph(ieph, 10);       % √a
e      = eph(ieph, 11);       % eccentricity

i_0    = eph(ieph, 12);       % inclination at ref (rad)
w      = eph(ieph, 13);       % argument of perigee (rad)
Omega_0= eph(ieph, 14);       % long. of asc. node (rad)
M_0    = eph(ieph, 15);       % mean anomaly (rad)

i_dot      = eph(ieph, 16);   % inclination rate (rad/s)
Omega_dot  = eph(ieph, 17);   % RAAN rate (rad/s)
delta_n    = eph(ieph, 18);   % mean motion diff. (rad/s)

C_uc = eph(ieph, 20);   C_us = eph(ieph, 21);
C_rc = eph(ieph, 22);   C_rs = eph(ieph, 23);
C_ic = eph(ieph, 24);   C_is = eph(ieph, 25);

% --- 2. 기본 항 ----------------------------------------------
a   = sqrt_a^2;                 % semi-major axis (m)
n_0 = sqrt(mu / a^3);           % computed mean motion
n   = n_0 + delta_n;            % corrected mean motion

t_k = gs - t_oe;                % 시간차
t_k = t_k - 604800*round(t_k/604800);   % ±302 400 s 범위로 wrap

M_k = M_0 + n*t_k;              % mean anomaly at tk
E_k = solveKeplerEq(M_k, e);    % eccentric anomaly (rad)
E_dot = n / (1 - e*cos(E_k));   % Ė

% --- 3. 보정 및 1계 미분 --------------------------------------
v_k  = atan2( sqrt(1-e^2)*sin(E_k), cos(E_k)-e );   % true anomaly
v_dot = sqrt(1-e^2)/(1 - e*cos(E_k))*E_dot;         % v̇   (ICD Eq.30)

phi_k  = v_k + w;                                   % argument of latitude
phi_dot = v_dot;

delta_u = C_us*cos(2*phi_k) + C_uc*sin(2*phi_k);
delta_r = C_rs*cos(2*phi_k) + C_rc*sin(2*phi_k);
delta_i = C_is*cos(2*phi_k) + C_ic*sin(2*phi_k);

u_dot = 2*(C_uc*cos(2*phi_k) - C_us*sin(2*phi_k))*phi_dot;
r_dot = a*e*sin(E_k)*E_dot + ...
        2*(C_rc*cos(2*phi_k) - C_rs*sin(2*phi_k))*phi_dot;
i_dot_corr = i_dot + ...
        2*(C_ic*cos(2*phi_k) - C_is*sin(2*phi_k))*phi_dot;

u_k = phi_k + delta_u;          % 보정 위도
u_dot = phi_dot + u_dot;        % u̇

r_k = a*(1 - e*cos(E_k)) + delta_r;   % 보정 반지름
i_k = i_0 + delta_i + i_dot*t_k;      % 보정 경사

% --- 4. 궤도면 좌표/속도 --------------------------------------
x_op  = r_k*cos(u_k);
y_op  = r_k*sin(u_k);

x_op_dot = r_dot*cos(u_k) - r_k*sin(u_k)*u_dot;
y_op_dot = r_dot*sin(u_k) + r_k*cos(u_k)*u_dot;

% --- 5. ECEF 회전 및 미분 -------------------------------------
Omega_k   = Omega_0 + (Omega_dot - we)*t_k - we*t_oe;
Omega_dot_corr = Omega_dot - we;

cO = cos(Omega_k);  sO = sin(Omega_k);
ci = cos(i_k);      si = sin(i_k);

vx =  x_op_dot*cO - y_op_dot*ci*sO ...
      - ( x_op*sO + y_op*ci*cO )*Omega_dot_corr;

vy =  x_op_dot*sO + y_op_dot*ci*cO ...
      + ( x_op*cO - y_op*ci*sO )*Omega_dot_corr;

vz =  y_op_dot*si + y_op*ci*i_dot_corr;

SatVel = [vx, vy, vz];          % (m/s)  |v| ≈ 3.9 km/s
end
