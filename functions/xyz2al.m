function al = xyz2al(rec_xyz, sat_xyz)
    % 위성 방향 벡터 (receiver -> satellite)
    vec = sat_xyz - rec_xyz;

    % 수신기 수직 방향 벡터 (Up)
    up = rec_xyz / norm(rec_xyz);

    % Elevation angle (단위: 도)
    el_rad = asin(dot(vec, up) / norm(vec));
    al = rad2deg(el_rad);
end