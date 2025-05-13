function jd = gwgs2jd(gw, gs)
    % GPS 에포크: 1980년 1월 6일 00:00:00 UTC -> JD = 2444244.5
    jd = 2444244.5 + 7 * gw + gs / 86400;
end
