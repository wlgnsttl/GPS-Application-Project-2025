function [dIno, lat_ipp, lon_ipp] = IonoGIM(Lat, Lon, TEC, vec_sat, vec_site, tc)

% 상수 정의
R_E = 6371000;
h_i = 450000;

% Z
vec_rho = vec_sat - vec_site;
Truellh = xyz2gd(vec_site);
[~, E] = xyz2azel(vec_rho, Truellh(1), Truellh(2));
Z = pi/2-E;
Z_p = asin(R_E * (sin(Z) / (R_E+h_i)));

% ipp
a = R_E;
b = R_E + h_i;
C = pi - (pi/2 + E + Z_p);

dis_ipp = sqrt(a^2 + b^2 - 2*a*b*cos(C));
vec_ipp = dis_ipp * vec_rho/norm(vec_rho);
vec_ipp = vec_ipp + vec_site;

llh_ipp = xyz2gd(vec_ipp);
lat_ipp = llh_ipp(1); lon_ipp = llh_ipp(2);

% interp
VTECgim = GetVTECgim(Lat, Lon, TEC, lat_ipp, lon_ipp, tc);

% VTEC -> STEC
OF = 1/cos(Z_p);
STEC = OF * VTECgim;

%dIno = 40.3 * STEC / (1575.42e6)^2;
dIno = STEC * 0.162372;