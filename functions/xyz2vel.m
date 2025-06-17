function TrueVel = xyz2vel(TruePos, dt_fine, dt_interp)
% XYZ2VEL - 위치 데이터를 dt_fine 간격으로 보간하고, dt_interp 간격으로 속도 복원
%
% 입력:
%   TruePos   (Nx4)     : [time, x, y, z] 형식의 위치 데이터 (원래는 1초 간격)
%   dt_fine   (scalar)  : 위치 보간 시간 간격 (ex. 0.1)
%   dt_interp (scalar)  : 속도 복원 시간 간격 (ex. 1.0)
%
% 출력:
%   TrueVel   (Mx3)     : dt_interp 간격의 속도 데이터 [vx, vy, vz]

    % 입력 시간 및 위치 분리
    time = TruePos(:,1);
    pos  = TruePos(:,2:4);

    % 보간 시간 벡터 생성
    time_fine = (time(1):dt_fine:time(end))';
    
    % 각 축 스플라인 보간
    xq = interp1(time, pos(:,1), time_fine, 'spline');
    yq = interp1(time, pos(:,2), time_fine, 'spline');
    zq = interp1(time, pos(:,3), time_fine, 'spline');

    % 속도 계산 (수치 미분)
    vx = gradient(xq, dt_fine);
    vy = gradient(yq, dt_fine);
    vz = gradient(zq, dt_fine);

    % 속도 복원 인덱스 계산
    step = round(dt_interp / dt_fine);
    idx  = 1:step:length(time_fine);

    % 속도 결과 추출
    TrueVel = [time_fine(idx), vx(idx), vy(idx), vz(idx)];
end
