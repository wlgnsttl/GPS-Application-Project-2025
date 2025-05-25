clear; close all;
clc;

addpath(genpath("data/"));
addpath(genpath("functions/"));

load('QM_RTAP5_250425_0748.mat');
load('eph_25115_1.mat');
% load('TruePos_ihvc.mat');
load('TruePos_RTAP5_250425_0748.mat');

TruePos = TruePos(:,2:4);
TrueVel = [0 0 0; diff(TruePos)];

[Lat,Lon,TEC] = ReadGIM('JPL0OPSFIN_20251150000_01D_02H_GIM.INX');

%% Date 정의

yyyy = 2025; 
mm = 3; 
dd = 20;

[gw, ~] = date2gwgs(yyyy, mm, dd, 0, 0, 0);
jd2mjd = 2400000.5;

%% 상수, 변수 정의

CCC = 299792458;
obsType = 103;

Truellh = xyz2gd(TruePos);

%% QM 선별

QM = SelectQM(arrQM, 100, obsType);
FinalTTs = unique(QM(:,1));

%% 추정에 필요한 반복 조건 및 초기값 설정

MaxIter = 5;
EpsStop = 1e-6; 
x = [TruePos(1,1:3), 1]; 
x = x';

%% 추정과정 시작

NoEpochs = length(FinalTTs);
estm = zeros(NoEpochs, 7);
nEst = 0;

for kE = 1:NoEpochs
    indxQM = find(QM(:,1) == FinalTTs(kE));
    QM1e = QM(indxQM, :);
    NoSats = size(QM1e, 1);
    gs = QM1e(1,1);

    for kIter = 1:MaxIter
        HTH = zeros(4,4);
        Hty = zeros(4,1);

        vec_site = x(1:3)';
        
        NoSatsUsed = 0;

        for kS = 1:NoSats
            prn = QM1e(kS, 2);         
            obs = QM1e(kS, 4);
            STT = obs/CCC;
            tc = gs - STT;
            ieph = PickEPH_multi(eph, prn, tc);
            toe = eph(ieph, 1); a = eph(ieph, 3); b = eph(ieph, 4); c = eph(ieph, 5); Tgd = eph(ieph, 6);
            
            [vec_sat, dRel] = getSatPos_lab(eph, ieph, tc);
            
            % 신호전달시간 보정
            vec_sat = RotSatPos(vec_sat, STT);

            vec_rho = vec_sat - vec_site;
            rho = norm(vec_rho);

            [az, el] = xyz2azel(vec_rho, Truellh(1), Truellh(2));

            % 고도각 7도 이하는 추정 과정에서 제거
            if el*180/pi <= 7
                continue;
            end
            
            % SV health이 0이 아니면 추정 과정에서 제거
            if eph(ieph, 19) > 0
                continue;
            end            

            % dIno = 0;
            [dIno, ipp_lat, ipp_lon] = IonoGIM(Lat, Lon, TEC, vec_sat, vec_site, tc);

            % dTrop = 0;
            dmjd = gwgs2jd(gw, tc) - jd2mjd;
            dlat = deg2rad(Truellh(1)); dlon = deg2rad(Truellh(2)); dhgt = Truellh(3);
            zd = pi/2 - el;
            [pres , temp, undu] = gpt_v1(dmjd, dlat, dlon, dhgt);
            % fprintf('%d\n', pres);
            [gmfh, ~] = gmf_f_hu(dmjd, dlat, dlon, dhgt, zd);
            dTrop = gmfh*(0.0022767*pres)/(1-0.00266*cos(2*dlat) - 0.00028*dhgt/1e3);

            % 보정값 합산
            dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel;
            com = rho + x(4) - CCC * dtSat + dIno + dTrop;
            y = obs - com;

            H  = [-vec_rho(1)/rho, -vec_rho(2)/rho, -vec_rho(3)/rho, 1];
            HTH = HTH + H'*H;
            Hty = Hty + H'*y;
            NoSatsUsed = NoSatsUsed + 1;
        end

        xhat = HTH \ Hty;
        x = x + xhat;

        if norm(xhat) < EpsStop
            nEst = nEst + 1;
            estm(nEst, 1) = gs;
            estm(nEst, 2:5) = x(1:4);
            estm(nEst, 6) = NoSats;
            estm(nEst, 7) = NoSatsUsed;
            vec_est = x(1:3)';
            break;
        end
    end
end
        
%% RMSE Calc
TTs = estm(:, 1);
XYZ = estm(:, 2:4);
NEV = xyz2topo2(XYZ, TruePos);
llh = xyz2gd(XYZ);

[rmse, horErr, ~, dim3Err] = nev2rmse(NEV);

%% Figure

close all; clc;
PlotPosRMSE(TTs, NEV, estm(:,6), estm(:,7));
figure;
geoplot(llh(:,1),llh(:,2))
%% Console disp
fprintf('%-15s : %6.3f [m]\n', 'Horizontal RMSE', rmse(1));
fprintf('%-15s : %6.3f [m]\n', 'Vertical RMSE', rmse(2));
fprintf('%-15s : %6.3f [m]\n\n', '3D RMSE', rmse(3));

fprintf('%-15s : %6.3f [m]\n', 'Max 2D Error', max(horErr));
fprintf('%-15s : %6.3f [m]\n', 'Max 3D Error', max(dim3Err));