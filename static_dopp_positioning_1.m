clear; close all;
clc;

addpath(genpath("data/"));
addpath(genpath("functions/"));

% load('QM_ihub_25079t.mat');
% load('eph_25079_1.mat');
% load('TruePos_ihub.mat');

% load('QM_ihvc_25079t.mat');
% load('eph_25079_1.mat');
load('TruePos_ihvc.mat');

load('QM_0529.mat');
eph = ReadEPH_multi('0529_1628250529_0728.nav');
%% 상수, 변수 정의

CCC = 299792458;
L1_lamda = 0.19029;

sys = 100;
obsType = 103;

Truellh = xyz2gd(TruePos);

%% QM 선별

QM = SelectQM(arrQM, sys, obsType);
FinalTTs = unique(QM(:,1));

%% 추정에 필요한 반복 조건 및 초기값 설정

MaxIter = 5;
EpsStop = 1e-5; 
x = [TruePos 0 0 0 0]; 
x = x';

%% 추정과정 시작

NoEpochs = length(FinalTTs);
estm = zeros(NoEpochs, 9);
nEst = 0;
MaxSnr = max(QM(:,7));

for kE = 1:NoEpochs
    idx   = QM(:,1)==FinalTTs(kE);
    QM1e  = QM(idx ,:);
    NoSats= size(QM1e,1);
    gs    = QM1e(1,1);
    
    for kIter = 1:MaxIter
        H = zeros(NoSats, 7);
        y = zeros(NoSats, 1);
    
        vec_rec_p = x(1:3);
        vec_rec_v = x(4:6);
        
        NoSatsUsed = 0;
        for kS = 1:NoSats
            prn = QM1e(kS,2);
            obs_dopp = QM1e(kS,6) * -L1_lamda;
    
            ieph  = PickEPH_multi(eph,prn,gs);
            
            if eph(ieph, 19) > 0
                continue;
            end  

            b = eph(ieph, 4);
            
            [vec_sat_p, ~] = getSatPos_lab(eph, ieph, gs);
            vec_sat_p = vec_sat_p';

            vec_sat_v = getSatVel(eph, ieph, gs)';
    
            vec_rho_p = vec_sat_p - vec_rec_p;
            rho = norm(vec_rho_p);
            h = vec_rho_p./rho;
    
            vec_rho_v = vec_sat_v - vec_rec_v;
            dr = h' * vec_rho_v;
            com = dr - CCC * b + x(7);

            g = -(1/ rho) * (eye(3) - h * h') * vec_rho_v;
            k = -h;

            NoSatsUsed = NoSatsUsed + 1;
            
            % 고도각 계산 및 저장
            [~, el] = xyz2azel(vec_rho_p, Truellh(1), Truellh(2));
            matrix_el(NoSatsUsed,:) = rad2deg(el);

            % snr 저장
            matrix_snr(NoSatsUsed,:) = QM1e(kS,7);

            H(NoSatsUsed,:) = [g', k', 1];
            y(NoSatsUsed) = obs_dopp - com;
        end
        
        if NoSatsUsed < 7
            continue;
        end

        H = H(1:NoSatsUsed, :);
        y = y(1:NoSatsUsed, :);
        
        % 고도각 가중치 
        W_el = 1;
        W_el = WeightEl(matrix_el);
        
        % SNR 가중치
        W_snr = 1;
        W_snr = WeightSNR(matrix_snr,MaxSnr);
        
        % 가충치 합산
        W = W_el .* W_snr;

        % 초기화
        matrix_el = 0;
        matrix_snr = 0;

        % xhat 계산
        xhat = pinv(H'*W*H)*H'*W*y;
        % xhat = pinv(H) * y;
        
        % 업데이트
        x = x + xhat;
    
        if norm(xhat) < EpsStop
            % fprintf("[%d] Err : %3.1fm\n", gs, norm(x(1:3) - TruePos(:)));

            nEst = nEst + 1;
            estm(nEst,1)   = gs;
            estm(nEst,2:8) = x;
            estm(nEst, 9) = NoSats;
            estm(nEst, 10) = NoSatsUsed;
            break
        end
    end
end

        
%% RMSE Calc

estm = estm(1:nEst, :);
estm_v = zeros(length(estm(:,1)),4);
estm_v(1,1:4) = estm(1,1:4);
estm_v(1,2:4) = TruePos(1,1:3);
for i = 1 : length(estm(:,1))-1
    estm_v(i+1,1) = estm(i+1,1);
    estm_v(i+1,2:4) = estm_v(i,2:4) + estm(i,5:7);
end

TTs = estm(:, 1);
XYZ = estm(:, 2:4);
VXYZ = estm(:, 5:7);
XYZ_v = estm_v(:, 2:4);

NEV = xyz2topo2(XYZ, TruePos);
NEV_v = xyz2topo2(XYZ_v, TruePos);

VNEV = xyz2topo2(VXYZ, [0 0 0]);

[rmse, horErr, verErr, dim3Err] = nev2rmse(NEV);
[rmse_v, horErr_v, verErr_v, dim3Err_v] = nev2rmse(NEV_v);

%% Figure

close all; clc;
PlotPosRMSE(TTs, NEV, estm(:,9), estm(:,10));
PlotVelRMSE(TTs, VNEV, estm(:,9), estm(:,10));
PlotPosRMSE(TTs, NEV_v, estm(:,9), estm(:,10));

%% Console disp

fprintf('%-15s : %6.3f [m]\n', 'Horizontal RMSE', rmse(1));
fprintf('%-15s : %6.3f [m]\n', 'Vertical RMSE', rmse(2));
fprintf('%-15s : %6.3f [m]\n\n', '3D RMSE', rmse(3));

fprintf('%-15s : %6.3f [m]\n', 'Max 2D Error', max(horErr));
fprintf('%-15s : %6.3f [m]\n', 'Max 3D Error', max(dim3Err));