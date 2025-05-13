function PlotSatVel(arrQM, eph, prn, obsType)
    % 지구 자전 각속도 벡터 [rad/s]
    omega_e = [0; 0; 7.2921159e-5];

    arrQM = SelectQM(arrQM, 0, obsType);
    arrQM = arrQM(arrQM(:,2) == prn);

    % 해당 PRN 위성의 시각만 추출
    TTs = arrQM(:, 1);
    
    % 결과 저장용
    result = zeros(length(TTs), 2);
    result(:,1) = TTs;

    % 각 시각에 대해 위치 및 속도 계산
    for i = 1:length(TTs)
        gs = TTs(i);
        ieph = PickEPH_multi(eph, prn, gs);  % ephemeris 선택

        % 위성 위치/속도 (ECEF 기준)
        r_ecef = getSatPos_lab(eph, ieph, gs);
        v_ecef = getSatVel(eph, ieph, gs);

        % 지구 자전 보정 → ECI 기준 속도
        v_eci = v_ecef(:) + cross(omega_e, r_ecef(:));
        result(i,2) = norm(v_eci);  % 공전 속도 크기
    end


    t = gs2gdhms(TTs);
    t = dms2deg(t(:,2),t(:,3),t(:,4));

    % 플로팅
    figure; hold on; grid on;
    plot(t, result(:,2), 'o-');
    xlim([0 24]);
    xlabel('hours');
    ylabel('ECI Velocity Magnitude [m/s]');
    title(sprintf('PRN [%d] ObsType [%d]', prn, obsType));
end
