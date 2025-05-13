function [nmea_data] = ReadNMEA(filename, option)
%GGA Data: 
%[UTC, latitude 위도, longitude 경도, gps_quality_indicator, nums_of_sat, HDOP, altitude, geoidal_separation, age_of_diff, diff_station_id]

    % 파일 열기
    fid = fopen(filename, 'r');
    
    if fid == -1
        error('파일을 열 수 없습니다.');
    end

    epoch_count = 0;
    nmea_data = [];

    while ~feof(fid)
        line = fgetl(fid);
        if ~isempty(line)
            if contains(line, 'GGA') && strcmp(option, 'GGA') %GGA 데이터만 추출
                epoch_count = epoch_count +1;
                GGA_data = readGGA(line);
                nmea_data(epoch_count, :) = GGA_data;
            elseif contains(line, 'GSA') && strcmp(option, 'GSA') %GGA 데이터만 추출
                epoch_count = epoch_count +1;
                GSA_data = readGSA(line);
                nmea_data(epoch_count, :) = GSA_data;
            end
        end
    end
    % 파일 닫기
    fclose(fid);
end

function [GGA_data] = readGGA(GGA_line)
%gga_line 슬라이싱
GGA_line = GGA_line(strfind(GGA_line, 'GGA'):end);
GGA_line = split(GGA_line, ',')';

%데이터 분류
HH = str2double(GGA_line{2}(1:2));
MM = str2double(GGA_line{2}(3:4));
SS = str2double(GGA_line{2}(5:end));
UTC = HH*3600 + MM*60 + SS;

%위도
latitude = str2double(GGA_line{3}(1:2)) + str2double(GGA_line{3}(3:end))/60;
if GGA_line{4} ~= 'N'
    latitude = -latitude;
end

%경도
longitude = str2double(GGA_line{5}(1:3)) + str2double(GGA_line{5}(4:end))/60;
if GGA_line{6} ~= 'E'
    longitude = -longitude;
end

gps_quality_indicator = str2double(GGA_line{7});
nums_of_sat = str2double(GGA_line{8});
HDOP = str2double(GGA_line{9});

%고도
altitude = str2double(GGA_line{10});
if strcmp(GGA_line{11}, 'cm')
    altitude = altitude/100;
end

%지오이드고
geoidal_separation = str2double(GGA_line{12});
if strcmp(GGA_line{13}, 'cm')
    geoidal_separation = geoidal_separation/100;
end

age_of_diff = str2double(GGA_line{14});
diff_station_id = str2double(extractBefore(GGA_line{15}, '*'));

GGA_data = [UTC, latitude, longitude, gps_quality_indicator, nums_of_sat, HDOP, altitude, geoidal_separation, age_of_diff, diff_station_id];
end

function [GSA_data] = readGSA(GSA_line)
GSA_line = GSA_line(strfind(GSA_line, 'GSA'):end);
GSA_line = split(GSA_line, ',')';

if GSA_line{2} == 'A'
    mode1 = 1;
else
    mode1 = 0;
end

mode2 = str2double(GSA_line{3});

prn = zeros(1, 12);
for i=1:12
    prn(i) = str2double(GSA_line{i+3});
end

VDOP = str2double(GSA_line{16});
HDOP = str2double(GSA_line{17});
PDOP = str2double(GSA_line{18});

GSA_data = [mode1, mode2, prn, VDOP, HDOP, PDOP];
end