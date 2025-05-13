function [dms] = deg2dms(deg)
    dd = fix(deg);
    mm = fix((deg-dd)*60);
    ss = ((deg-dd)*60 - mm)*60;
    [dms] = [dd, mm, ss];