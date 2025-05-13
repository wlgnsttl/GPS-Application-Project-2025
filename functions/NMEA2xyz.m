function [xyz] = NMEA2xyz(filepath)

GGA = ReadNMEA(filepath, "GGA");

TTs = GGA(:,1);

LatLon = GGA(:,2:3);
H = GGA(:,7) + GGA(:,8);

LLH = [LatLon, H];

xyz = [TTs, gd2xyz(LLH)];