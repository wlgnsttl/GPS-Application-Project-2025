function [gps_week_number, gps_week_day, gps_week_seconds] = date2gwgdgs(yyyy,mm,dd,HH,MM,SS)
    %date(yyyy-mm-dd HH:MM:SS)를 입력받아
    %GPS Week Number, GPS Week Day, GPS Week Seconds 계산
    
    %date를 JD로 변환
    jd = date2jd(yyyy,mm,dd,HH,MM,SS);
    
    %GPS Week Number계산
    gps_week_number = floor((jd - 2444244.5)/7);
    
    %N은 요일(N=0은 Monday, N=1은 화요일)
    N = mod(floor(jd + 0.5), 7);

    %N을 GPS Week Day로 변환(gps_week_day = 0은 일요일)
    gps_week_day = mod(N + 1, 7);
    
    %GPS Week Seconds 계산
    gps_week_seconds = 86400*gps_week_day + HH*3600 + MM*60 + SS;
end