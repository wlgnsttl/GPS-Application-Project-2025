clear; close all;
clc;
%%
addpath(genpath("data/"));
addpath(genpath("functions/"));

load('QM_RTAP5_250425_0748.mat');
load('eph_25115_1.mat');
load('TruePos_RTAP5_250425_0748.mat');

[Lat,Lon,TEC] = ReadGIM('JPL0OPSFIN_20251150000_01D_02H_GIM.INX');

obsType = 103;

QM = SelectQM(arrQM, 100, obsType);
FinalTTs = unique(QM(:,1));

TruePos = TruePos(:,2:4);

%% 상수
CCC = 299792458; % 광속 (m/s)
Truellh = xyz2gd(TruePos);
%% Date 정의

yyyy = 2025; 
mm = 4; 
dd = 25;

[gw, ~] = date2gwgs(yyyy, mm, dd, 0, 0, 0);
jd2mjd = 2400000.5;
%% 초기 상태 및 공분산
x_hat = [TruePos(1,1:3)';0; 0; 0; 0]; % [x; y; z; b] 초기 위치 추정
P = blkdiag(eye(3),eye(3),100) ;     % 초기 공분산 행렬
% Q = diag([0, 0, 0, 100]);
Q = blkdiag(eye(3),eye(3),1e2);
% Q = Q.*Q;
A = eye(7);
%% 측정 노이즈
estm =[];
for epoch = 1:length(FinalTTs)
    gs = FinalTTs(epoch);
    % 현재 epoch에서 사용 가능한 관측값
    idx = QM(:,1) == gs;
    QM_obs = QM(idx, :);
    NoSats = QM_obs(:,2);
    
    % 예측 단계
    xp = A*x_hat;
    Pp = A*P*A' + Q;
    
    % 측정 예측 및 칼만 필터 업데이트
    % H = zeros(length(NoSats), 4);
    % h = zeros(length(NoSats), 1);
    % R = eye(length(NoSats));
    clear H com R
    obs = [];
    R_el = [];
    NosatsUsed=1;
    for i = 1:length(NoSats)
        obs_t = QM_obs(i,4); % pseudorange [m]
        STT = obs_t/CCC;
        prn = NoSats(i);
        tc = QM_obs(i,1) - STT;
        ieph  = PickEPH_multi(eph,prn,gs);
        toe = eph(ieph, 1); a = eph(ieph, 3); b = eph(ieph, 4); c = eph(ieph, 5); Tgd = eph(ieph, 6);

        [vec_sat, dRel] = getSatPos_lab(eph, ieph, tc);
        vec_sat = RotSatPos(vec_sat,STT);

        vec_rho = vec_sat' - xp(1:3,1);  % 위성 - 수신기
        rho = norm(vec_rho);

        [az, el] = xyz2azel(vec_rho', Truellh(1),Truellh(2));
        if rad2deg(el)<7
            continue
        end
        if eph(ieph, 19) > 0
                continue;
        end
        % dIno = 0;
        [dIno, ipp_lat, ipp_lon] = IonoGIM(Lat, Lon, TEC, vec_sat, xp(1:3)', tc);

        % dTrop = 0;
        dmjd = gwgs2jd(gw, tc) - jd2mjd;
        dlat = deg2rad(Truellh(1)); dlon = deg2rad(Truellh(2)); dhgt = Truellh(3);
        zd = pi/2 - el;
        [pres , temp, undu] = gpt_v1(dmjd, dlat, dlon, dhgt);

        [gmfh, ~] = gmf_f_hu(dmjd, dlat, dlon, dhgt, zd);
        dTrop = gmfh*(0.0022767*pres)/(1-0.00266*cos(2*dlat) - 0.00028*dhgt/1e3);
        % 보정값 합산
        dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel;

               
        obs(NosatsUsed,1) = obs_t;
        com(NosatsUsed,1) = rho + xp(4) - CCC * dtSat + dIno + dTrop;
        H(NosatsUsed,:) = [-vec_rho(1)/rho, -vec_rho(2)/rho, -vec_rho(3)/rho,1];
        R_el(NosatsUsed,1) =(el);
        NosatsUsed= NosatsUsed+1;
    end
    for i = 1 :length(NoSats)
        obs_t = QM_obs(i,4); % pseudorange [m]
        STT = obs_t/CCC;
        prn = NoSats(i);
        tc = QM_obs(i,1) - STT;
        ieph  = PickEPH_multi(eph,prn,gs);
        toe = eph(ieph, 1); a = eph(ieph, 3); b = eph(ieph, 4); c = eph(ieph, 5); Tgd = eph(ieph, 6);

        [vec_sat_p, vec_sat_v] = getSatVel_diff(eph, ieph, tc, STT);

        vec_rho_p = vec_sat_p - xp(1:3,1);
        rho = norm(vec_rho_p);
        h = vec_rho_p./rho;
        [az, el] = xyz2azel(vec_rho_p', Truellh(1),Truellh(2));
        if rad2deg(el)<7
            continue
        end
        if eph(ieph, 19) > 0
                continue;
        end


        vec_rho_v = vec_sat_v - xp(4:6,1);
        dr = h' * vec_rho_v;
        com = dr - CCC * b + xp(7);

        % g = -(1/ rho) * (eye(3) - h * h') * vec_rho_v;
        k = -h;
        
    end
    
    % 칼만 이득
    % R = diag(1 ./ sin(R_el).^2); 
    % R = eye(k-1).*3;
    R = WeightEl(1./R_el);
    y = obs - com;
    K = Pp * H' / (H * Pp * H' + R);

    % 상태 갱신
    x_hat = xp + K * (y);
    P = Pp - K*H*Pp;

    
    % 결과 저장
    estm = [estm; gs, x_hat',NosatsUsed-1,length(NoSats)];
end

%% RMSE Calc
idx = estm(:, 1) ~= 0;
estm = estm(idx, :); TruePos = TruePos(idx, :);
%% Figure
TTs = estm(:,1);
NEV = xyz2topo2(estm(:,2:4),TruePos);
close all; clc;
PlotPosRMSE(TTs, NEV, estm(:,7), estm(:,6));