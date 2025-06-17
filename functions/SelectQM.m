function QM = SelectQM(arrQM, sys, obsType)
    QM = arrQM;  % 초기화

    if sys ~= 0
        QM = QM(any(floor(QM(:,2) / 100) == (sys/100), 2), :);
    end

    if obsType ~= 0
        QM = QM(any(QM(:,3) == obsType, 2), :);
    end
end
