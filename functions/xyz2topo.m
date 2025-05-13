function [dnev] = xyz2topo(dxyz, lat, longi)
%lat,longi(deg)에서의 dxyz를 입력받아 dnev로 반환

%deg2rad
lat = lat*pi/180;
longi = longi*pi/180;

R = [-sin(lat)*cos(longi), -sin(lat)*sin(longi), cos(lat);
    -sin(longi), cos(longi), 0;
    cos(lat)*cos(longi), cos(lat)*sin(longi), sin(lat)];

dnev = (R*dxyz')';