function [dNEV] = xyz2topo(dXYZ, lat_deg, lon_deg)
% ECEF 차이 벡터를 NEV 좌표계로 변환
% 입력:
%   dXYZ     : nx3 차이 벡터
%   lat_deg  : 기준 위도 [deg]
%   lon_deg  : 기준 경도 [deg]
% 출력:
%   dNEV     : nx3 NEV 벡터 (North, East, Vertical)

    % 도 → 라디안
    lat = deg2rad(lat_deg);
    lon = deg2rad(lon_deg);

    % NEV 회전 행렬 (ECEF → NEV)
    R = [-sin(lat)*cos(lon), -sin(lat)*sin(lon),  cos(lat);
         -sin(lon),           cos(lon),           0;
          cos(lat)*cos(lon),  cos(lat)*sin(lon),  sin(lat)];

    % 회전 적용
    dNEV = (R * dXYZ')';
end
