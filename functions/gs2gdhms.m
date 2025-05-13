function [dhms, H] = gs2gdhms(gs)
%gs를 [GPS Week Day H M S] 로 변환

gs_length = length(gs);
dhms = NaN(gs_length, 4);
day_second = 86400;

for i=1:gs_length
    temp = gs(i);
    gd = floor(temp/day_second);
    temp = temp - gd*day_second;

    h = floor(temp/3600);
    temp = temp-h*3600;

    m = floor(temp/60);
    s = temp-m*60;

    dhms(i,:) = [gd, h, m, s];
end

H = dms2deg(dhms(:,2), dhms(:,3), dhms(:,4));