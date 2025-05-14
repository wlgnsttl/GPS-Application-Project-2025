function [dNEV] = xyz2topo3(XYZ, TruePos)
% [dNEV] = xyz2topo3(XYZ, truePos) 
% 이동측위 version
% 
% dNEV    : truePos를 원점으로 계산된 nx3 NEV 좌표 행렬
% XYZ     : nx3의 XYZ(ECEF) 행렬
% truePos : dNEV의 원점이 되는 행벡터 xyz 좌표(Matrix)

trueLLH = xyz2gd(TruePos);
dXYZ = XYZ - TruePos;

nEph = size(dXYZ, 1);
dNEV = zeros(nEph, 3);

for iter = 1:nEph
    dnev = xyz2topo(dXYZ(iter, :), trueLLH(iter,1), trueLLH(iter,2));
    dNEV(iter,:) = dnev;
end

