function [xyz] = NMEA2xyz(filepath)

GGA = ReadNMEA(filepath, 'GGA');

TTs = GGA(:,1);

LL = GGA(:, 2:3);
H = GGA(:,7) + GGA(:,8);

LLH = [LL, H];

pos = gd2xyz(LLH);

xyz = [TTs pos];