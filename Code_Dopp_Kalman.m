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
TrueVel = [0 0 0; diff(TruePos)];

%% 상수
CCC = 299792458; % 광속 (m/s)
L1_lamda = 0.19029;
Truellh = xyz2gd(TruePos);
%% Date 정의

yyyy = 2025; 
mm = 4; 
dd = 25;

[gw, ~] = date2gwgs(yyyy, mm, dd, 0, 0, 0);
jd2mjd = 2400000.5; 

% 초기 상태 및 공분산
x_hat = [TruePos(1,1:3)'; 0; 0; 0; 0; 0]; % [x; y; z; b] 초기 위치 추정
P = blkdiag(eye(3).*0.1,10000,eye(3),1000);     % 초기 공분산 행렬
% P(1,5) = 1;
% P(2,6) = 1;
% P(3,7) = 1;
% P(5,1) = 1;
% P(6,2) = 1;
% P(7,3) = 1;
% Q = diag([0, 0, 0, 100]);
Q = blkdiag(eye(3)*0.3,3e+2,eye(3)*17,9e+1);
Q = Q.*Q;
%x,y,z 1m~5m
%vx, vy, vz 0.1m 0.5m
A = eye(8);
A(1,5) = 1;
A(2,6) = 1;
A(3,7) = 1;
%% 측정 노이즈
estm =[];
for epoch = 1:length(FinalTTs)
    gs = FinalTTs(epoch);
    % 현재 epoch에서 사용 가능한 관측값
    idx = QM(:,1) == gs;
    QM_obs = QM(idx, :);
    NoSats = QM_obs(:,2);
    MaxSnr = max(QM_obs(:,7));
    
    % 예측 단계
    xp = A*x_hat;
    Pp = A*P*A' + Q;
    
    % 측정 예측 및 칼만 필터 업데이트
    % H = zeros(length(NoSats), 4);
    % h = zeros(length(NoSats), 1);
    % R = eye(length(NoSats));
    clear H com R obs_c obs_d com_c com_d H_c H_d
    obs = [];
    R_el = [];
    R_snr = [];
    NoSatsUsed = 0;
    for i = 1:length(NoSats)
        obs_t = QM_obs(i,4); % pseudorange [m]
        obs_dopp = QM_obs(i,6) * -L1_lamda;
        STT = obs_t/CCC;
        prn = NoSats(i);
        tc = QM_obs(i,1) - STT;
        snr =  QM_obs(i,7);
        ieph  = PickEPH_multi(eph,prn,gs);
        toe = eph(ieph, 1); a = eph(ieph, 3); b = eph(ieph, 4); c = eph(ieph, 5); Tgd = eph(ieph, 6);
        % pseudorange
        [vec_sat, dRel] = getSatPos_lab(eph, ieph, tc);
        vec_sat = RotSatPos(vec_sat,STT);

        vec_rho = vec_sat' - xp(1:3,1);  % 위성 - 수신기
        rho = norm(vec_rho);

        [az, el] = xyz2azel(vec_rho', Truellh(epoch,1),Truellh(epoch,2));
        if rad2deg(el)<10
            continue;
        end
        if eph(ieph, 19) > 0
            continue;
        end
        if QM_obs(i,7) < 30
            continue;
        end
        % dIno = 0;
        [dIno, ipp_lat, ipp_lon] = IonoGIM(Lat, Lon, TEC, vec_sat, xp(1:3)', tc);

        % dTrop = 0;
        dmjd = gwgs2jd(gw, tc) - jd2mjd;
        dlat = deg2rad(Truellh(epoch,1)); dlon = deg2rad(Truellh(epoch,2)); dhgt = Truellh(epoch,3);
        zd = pi/2 - el;
        [pres , temp, undu] = gpt_v1(dmjd, dlat, dlon, dhgt);

        [gmfh, ~] = gmf_f_hu(dmjd, dlat, dlon, dhgt, zd);
        dTrop = gmfh*(0.0022767*pres)/(1-0.00266*cos(2*dlat) - 0.00028*dhgt/1e3);
        % 보정값 합산
        dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel;
        % doppler
        b = eph(ieph, 4);
        [vec_sat_p, vec_sat_v] = getSatVel_diff(eph, ieph, tc, STT);

        vec_rho_p = vec_sat_p - xp(1:3);
        rho = norm(vec_rho_p);
        h = vec_rho_p./rho;

        vec_rho_v = vec_sat_v - xp(5:7);
        dr = h' * vec_rho_v;
        
        g = -(1/ rho) * (eye(3) - h * h') * vec_rho_v;
        k_d = -h;

        NoSatsUsed = NoSatsUsed+1;
        obs_c(NoSatsUsed,1) = obs_t;
        com_c(NoSatsUsed,1) = rho + xp(4) - CCC * dtSat + dIno + dTrop;
        obs_d(NoSatsUsed,1) = obs_dopp;
        com_d(NoSatsUsed,1) = dr - CCC * b + xp(8);
        H_c(NoSatsUsed,:) = [-vec_rho(1)/rho, -vec_rho(2)/rho, -vec_rho(3)/rho, 1, 0, 0, 0, 0];
        H_d(NoSatsUsed,:) = [0, 0, 0, 0, -vec_rho(1)/rho, -vec_rho(2)/rho, -vec_rho(3)/rho, 1];
        R_el(NoSatsUsed,1) =el;
        R_snr(NoSatsUsed,1) =snr;
    end
    obs(1:NoSatsUsed,1) = obs_c;
    obs(NoSatsUsed+1:2*NoSatsUsed,1) = obs_d;
    com(1:NoSatsUsed,1) = com_c;
    com(NoSatsUsed+1:2*NoSatsUsed,1) = com_d;
    H(1:NoSatsUsed,:) = H_c;
    H(NoSatsUsed+1:2*NoSatsUsed,:) = H_d;
    % 칼만 이득
    R_el = [R_el*100; R_el];
    R_snr = [R_snr*100; R_snr];
    % R = WeightEl(1./R_el);% + WeightSNR(1./R_snr,MaxSnr);
    % R = eye(2*NoSatsUsed).*2^2;
    R = blkdiag(eye(NoSatsUsed)*2.5^2,eye(NoSatsUsed)*0.05^2);
    % R = blkdiag(eye(NoSatsUsed)*2.5^2);
    % R 값 3m^2으로 주기
    y = obs - com;
    K = Pp * H' / (H * Pp * H' + R);

    % 상태 갱신
    x_hat = xp + K * (y);
    P = Pp - K*H*Pp;

    
    % 결과 저장
    estm = [estm; gs, x_hat',NoSatsUsed,length(NoSats)];
end

%% Figure
TTs = estm(:,1);
XYZ = estm(:,2:4);
llh = xyz2gd(XYZ);
NEV = xyz2topo2(estm(:,2:4),TruePos);

close all; clc;
PlotPosRMSE(TTs, NEV, estm(:,11), estm(:,10));
figure;
hold on;
grid on;
scatter(NEV(:,1),NEV(:,2),'r.')
xlabel('N');
ylabel('E');
hold off;

figure;
geoplot(llh(:,1),llh(:,2),'bo-','MarkerSize',5);
hold on;
geoplot(Truellh(:,1), Truellh(:,2),'ro-','MarkerSize',5);