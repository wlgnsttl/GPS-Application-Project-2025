function deltaR = sagnac(rec_xyz, sat_xyz)
% Sagnac 효과 보정값 계산 (단위: 미터)
% 입력:
%   rec_xyz : 수신기 ECEF 좌표 [x; y; z] (단위: m)
%   sat_xyz : 위성 ECEF 좌표 [x; y; z] (단위: m)
% 출력:
%   deltaR : Sagnac 효과로 인한 거리 보정값 (단위: m)

    % 지구 자전 각속도 벡터 (rad/s)
    omega_e = 7.2921151467e-5; % rad/s
    Omega = [0; 0; omega_e];

    % 거리 벡터
    r_rs = sat_xyz(:) - rec_xyz(:); % 위성 - 수신기

    % 보정량 계산: deltaR = (Omega × sat_xyz) · r_rs / c
    cross_term = cross(Omega, sat_xyz(:));
    c = 299792458; % 빛의 속도 (m/s)

    deltaR = dot(cross_term, r_rs) / c;
end

