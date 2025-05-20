function [el_rad, az_rad] = xyz2elaz(rec_xyz, sat_xyz)
% 입력:
%   rec_xyz : 수신기 ECEF 좌표 (3x1)
%   sat_xyz : 위성 ECEF 좌표 (3x1)
% 출력:
%   el_rad  : 고도각 (radian)
%   az_rad  : 방위각 (radian)

% 거리 벡터
vec_rs = sat_xyz(:) - rec_xyz(:);

% 수신기 위치에서 ENU 기준축 정의
x = rec_xyz(1);
y = rec_xyz(2);
z = rec_xyz(3);

r = norm(rec_xyz);
phi = asin(z / r);     % geocentric latitude (radian)
lambda = atan2(y, x);  % longitude (radian)

% ENU 기준 벡터 (geocentric 기준)
east  = [-sin(lambda);  cos(lambda); 0];
north = [-sin(phi)*cos(lambda); -sin(phi)*sin(lambda); cos(phi)];
up    = [ cos(phi)*cos(lambda);  cos(phi)*sin(lambda); sin(phi)];

% 거리 벡터를 ENU 좌표로 투영
e = dot(vec_rs, east);
n = dot(vec_rs, north);
u = dot(vec_rs, up);

% 고도각
el_rad = atan2(u, sqrt(e^2 + n^2));

% 방위각 (북 기준, 시계방향)
az_rad = atan2(e, n);
az_rad = mod(az_rad, 2*pi); % 0 ~ 2pi 범위
end