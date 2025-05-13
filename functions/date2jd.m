function jd = date2jd(yyyy,mm,dd,HH,MM,SS)
    %date(yyyy-mm-dd HH:MM:SS)를 Julian Date(JD)로 변환합니다.

    if mm <= 2
        y = yyyy-1;
        m = mm+12;
    else
        y = yyyy;
        m = mm;
    end
    
    ut = HH + MM/60 + SS/3600; %UT계산
    jd = floor(365.25*y) + floor(30.6001*(m+1)) + dd + ut/24 + 1720981.5; %JD계산
end
