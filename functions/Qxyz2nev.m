function Q_nev = Qxyz2nev(Q_xyz4x4, lat_deg, lon_deg)
% Qxyz2nev - ECEF 4x4 공분산 행렬을 NEV 기준으로 변환
%
% INPUT:
%   Q_xyz4x4 : 4x4 공분산 행렬 (위치 3 + 클럭 바이어스 1)
%   lat_deg  : 수신기 위도 [deg]
%   lon_deg  : 수신기 경도 [deg]
%
% OUTPUT:
%   Q_nev : 4x4 공분산 행렬 (NE, E, V, clock bias)

% 위도/경도 (rad)
lat = deg2rad(lat_deg);
lon = deg2rad(lon_deg);

% 3x3 회전 행렬 (ECEF → NEV)
T33 = [-sin(lat)*cos(lon), -sin(lat)*sin(lon),  cos(lat);
       -sin(lon),           cos(lon),           0;
       -cos(lat)*cos(lon), -cos(lat)*sin(lon), -sin(lat)];

% 4x4 회전 행렬 (마지막 행/열은 clock bias → 그대로 유지)
T44 = eye(4);
T44(1:3,1:3) = T33;

% 공분산 행렬 변환
Q_nev = T44 * Q_xyz4x4 * T44';

end
