function [rmse, horErr, verErr, dim3Err] = nev2rmse(dNEV)
    % nev matrix(n by [n, e, v]) 
    % -> 
    % rmse [Horizontal RMSE, Vertical RMSE, 3D RMSE]
    % 
    % -------------Expand Option---------------
    % Horizontal Error(n by err)
    % Vertical Error(n by err)
    % 3D Error(n by err)

    % 수평 오차: 첫 두 열(x, y)의 유클리드 노름 계산
    horErr = sqrt(sum(dNEV(:,1:2).^2, 2));
    
    % 수직 오차: 세 번째 열(z)
    verErr = abs(dNEV(:,3));
    
    % 3차원 오차: 첫 세 열(x, y, z)의 유클리드 노름 계산
    dim3Err = sqrt(sum(dNEV(:,1:3).^2, 2));
    
    % 각각의 RMSE 계산
    rmse = [sqrt(mean(horErr.^2)), sqrt(mean(verErr.^2)), sqrt(mean(dim3Err.^2))];
end
