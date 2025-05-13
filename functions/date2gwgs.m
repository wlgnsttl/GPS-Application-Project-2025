function [gw, gs] = date2gwgs(year, month, day, hh, mm, ss)
    %Calendar Day와 시간(시, 분, 초)로 GPS Week Number와 GPS Week Seconds를 계산

    [gw, ~, gs] = date2gwgdgs(year, month, day, hh, mm, ss);
end

