function [dNEV] = xyz2topo2(XYZ, truePos)
% [dNEV] = xyz2topo2(XYZ, truePos)
% 
% dNEV    : truePos를 원점으로 계산된 nx3 NEV 좌표 행렬
% XYZ     : nx3의 XYZ(ECEF) 행렬
% truePos : dNEV의 원점이 되는 행벡터 xyz 좌표

trueLLH = xyz2gd(truePos);
dXYZ = XYZ - truePos;
dNEV = xyz2topo(dXYZ, trueLLH(1), trueLLH(2));