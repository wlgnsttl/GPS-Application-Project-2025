function [xyz] = gd2xyz(gd)
% 경위도 좌표를 WGS-84제원으로 

% WGS-84 제원
a = 6378137.0;    % 반지름
f = 1 / 298.257223563;  % 편평도
b = a * (1 - f);  % 치수
e = sqrt((a^2 - b^2) / a^2);  % 이심률

xyz = NaN(size(gd, 1), 3);

for i = 1:size(gd, 1)
    % 위도와 경도 값 확인
    lat = gd(i, 1);  % 위도
    longi = gd(i, 2);  % 경도
    h = gd(i, 3);    % 고도
    
    % 위도와 경도에 대해 계산
    N = (a^2) / sqrt(a^2 * cosd(lat)^2 + b^2 * sind(lat)^2);
    X = (N + h) * cosd(lat) * cosd(longi);
    Y = (N + h) * cosd(lat) * sind(longi);
    Z = (N * (1 - e^2) + h) * sind(lat);
    
    xyz(i, :) = [X, Y, Z];  % 결과 저장
end

end
