function QM = SelectQM(arrQM, sys, obsType)
    QM = arrQM;  % 초기화

    if sys ~= 0
        QM = QM(floor(QM(:,2) / 100) == (sys/100), :);
    end

    if obsType ~= 0
        QM = QM(QM(:,3) == obsType, :);
    end
end
