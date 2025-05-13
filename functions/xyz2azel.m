function [az, el] = xyz2azel(vec_rho, lat, lon)
% xyz2azel : 좌표 변환 함수
% vec_rho : ECEF 좌표계상의 벡터 (예: [x, y, z])
% lat, lon: 기준점의 위도와 경도 (단위: deg)
%
% 출력:
% az : 방위각 (azimuth, 라디안)
% el : 고도각 (elevation, 라디안)

% 위도와 경도 할당
phi = lat; 
lam = lon;

% 변환 행렬 A 구성 (ECEF -> ENU 또는 NEV 좌표계)
A = [ -sind(lam),              cosd(lam),              0;
     -sind(phi)*cosd(lam),     -sind(phi)*sind(lam),     cosd(phi);
      cosd(phi)*cosd(lam),      cosd(phi)*sind(lam),     sind(phi) ];

% 입력 벡터를 열 벡터로 변환 후 변환 행렬 적용
Y = A * vec_rho(:);

rho_e = Y(1);
rho_n = Y(2);
rho_v = Y(3);

% atan2를 사용하여 올바른 사분면 고려 (방위각 계산)
az = atan2(rho_e, rho_n);
az = mod(az, 2*pi);

% 고도각 계산: 수평면에 대한 각도, atan2 사용
el = atan2(rho_v, sqrt(rho_e^2 + rho_n^2));
