clear; close all;
clc;
%%
addpath(genpath("data/"));
addpath(genpath("functions/"));

load('TruePos_ihvc.mat');
load('QM_ihvc_25079t.mat');
load('eph_25079_1.mat');

[Lat,Lon,TEC] = ReadGIM('JPL0OPSFIN_20250790000_01D_02H_GIM.INX');

obsType = 103;

QM = SelectQM(arrQM, 100, obsType);
FinalTTs = unique(QM(:,1));



%% 상수
CCC = 299792458; % 광속 (m/s)
Truellh = xyz2gd(TruePos);
%% Date 정의

yyyy = 2025; 
mm = 3; 
dd = 20;

[gw, ~] = date2gwgs(yyyy, mm, dd, 0, 0, 0);
jd2mjd = 2400000.5;
%% 초기 상태 및 공분산
x_hat = [TruePos'; 0]; % [x; y; z; b] 초기 위치 추정
P = blkdiag(eye(3),100) ;     % 초기 공분산 행렬
% Q = diag([0, 0, 0, 100]);
Q = blkdiag(eye(3)*1e-6,1e2).*100000;
% Q = Q.*Q;
A = eye(4);
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
    k=1;
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
        
           
        obs(k,1) = obs_t;
        com(k,1) = rho + xp(4) - CCC * dtSat + dIno + dTrop;
        H(k,:) = [-vec_rho(1)/rho, -vec_rho(2)/rho, -vec_rho(3)/rho, 1];
        R_el(k,1) =(el);
        k= k+1;
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
    estm = [estm; gs, x_hat',k-1,length(NoSats)];
end

%% Figure
TTs = estm(:,1);
NEV = xyz2topo2(estm(:,2:4),TruePos);
close all; clc;
PlotPosRMSE(TTs, NEV, estm(:,7), estm(:,6));