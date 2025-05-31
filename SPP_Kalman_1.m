clear; close all;
clc;
%%
addpath(genpath("data/"));
addpath(genpath("functions/"));
load('TruePos_ihvc.mat');

load('QM_ihvc_25079t.mat');
load('eph_25079_1.mat');
obsType = 103;

QM = SelectQM(arrQM, 100, obsType);
FinalTTs = unique(QM(:,1));



%% 상수
CCC = 299792458; % 광속 (m/s)
Truellh = xyz2gd(TruePos);
%% 초기 상태 및 공분산
x_hat = [TruePos'; 0]; % [x; y; z; b] 초기 위치 추정
P = blkdiag(eye(3)*1,100) ;     % 초기 공분산 행렬
% Q = diag([0, 0, 0, 100]);
Q = blkdiag(eye(3)*1e-6,1e2);
% Q = Q.*Q;
A = eye(4);
%% 측정 노이즈
estm = zeros(1,6);

for epoch = 1:length(FinalTTs)
    t = FinalTTs(epoch);
    % 현재 epoch에서 사용 가능한 관측값
    idx = QM(:,1) == t;
    QM_obs = QM(idx, :);
    NoSats = QM_obs(:,2);
    
    % 예측 단계
    xp = A*x_hat;
    Pp = A*P*A' + Q;
    
    % 측정 예측 및 칼만 필터 업데이트
    % H = zeros(length(NoSats), 4);
    % h = zeros(length(NoSats), 1);
    % R = eye(length(NoSats));
    clear H h R
    obs = [];
    R_el = [];
    k=1;
    for i = 1:length(NoSats)
        obs_t = QM_obs(i,4); % pseudorange [m]
        STT = obs_t/CCC;
        prn = NoSats(i);
        tc = QM_obs(i,1) - STT;
        ieph  = PickEPH_multi(eph,prn,t);

        [sat_pos, dRel] = getSatPos_lab(eph, ieph, tc);
        sat_pos = RotSatPos(sat_pos,STT);

        vec_rho = sat_pos' - xp(1:3,1);  % 위성 - 수신기
        rho = norm(vec_rho);

        [~, el] = xyz2azel(-vec_rho', Truellh(1),Truellh(2));
        if rad2deg(el)<7
            continue
        end
        obs(k,1) = obs_t;
        h(k,1) = rho + xp(4) - CCC * dRel;
        H(k,:) = [-vec_rho'/rho, 1]; 
        R_el(k,1) =(el);
        k= k+1;
    end
    
    % 칼만 이득
    % R = diag(1 ./ sin(R_el).^2); 
    R = eye(k-1).*3;
    y_resid = obs - h;
    K = Pp * H' / (H * Pp * H' + R);
    
    % 상태 갱신
    x_hat = xp + K * (y_resid);
    P = Pp - K*H*Pp;
    
    % 결과 저장
    estm = [estm; t, x_hat',length(NoSats)];
end

%% Figure
TTs = estm(:,1);
NEV = xyz2topo2(estm(:,2:4),TruePos);
close all; clc;
PlotPosRMSE(TTs, NEV, estm(:,6), estm(:,6));