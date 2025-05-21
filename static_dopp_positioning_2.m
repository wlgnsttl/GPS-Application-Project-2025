clear; close all;
clc;

addpath(genpath("data\"));
addpath(genpath("functions\"));

load('QM_ihub_25079t.mat');
load('eph_25079_1.mat');
load('TruePos_ihub.mat');

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
x = [TruePos 0]; 
x = x';

%% 추정과정 시작

NoEpochs = length(FinalTTs);
estm = zeros(NoEpochs, 7);
nEst = 0;

for kE = 1:NoEpochs
    idx   = QM(:,1)==FinalTTs(kE);
    QM1e  = QM(idx ,:);
    NoSats= size(QM1e,1);
    gs    = QM1e(1,1);
    
    for kIter = 1:MaxIter
        H = zeros(NoSats, 4);
        y = zeros(NoSats, 1);
    
        vec_rec_p = x(1:3);
        vec_rec_v = [0 0 0]';
        
        NoSatsUsed = 0;
        for kS = 1:NoSats
            prn = QM1e(kS,2);
            obs_dopp = QM1e(kS,6) * -L1_lamda;
    
            ieph  = PickEPH_multi(eph,prn,gs);
            
            if eph(ieph, 19) > 0
                continue;
            end  

            b = eph(ieph, 4); % 위성 시계오차 변화율
            
            [vec_sat_p, ~] = getSatPos_lab(eph, ieph, gs);
            vec_sat_p = vec_sat_p';

            vec_sat_v = getSatVel(eph, ieph, gs)'; % Get Sat Velocity
    
            vec_rho_p = vec_rec_p - vec_sat_p;
            rho = norm(vec_rho_p);
            h = vec_rho_p./rho;
 
        
            vec_rho_v = vec_rec_v - vec_sat_v;
            dr = h' * vec_rho_v;
            com = dr - CCC * b + x(4);

            gt = (vec_rho_v' / rho) * (eye(3) - h * h');

            NoSatsUsed = NoSatsUsed + 1;
            
            % 고도각 계산 및 저장
            [~, el] = xyz2azel(vec_rho_p, Truellh(1), Truellh(2));
            matrix_el(NoSatsUsed,:) = rad2deg(el);
            
            H(NoSatsUsed,:) = [gt, 1];
            y(NoSatsUsed) = obs_dopp - com;
        end
        
        if NoSatsUsed < 4
            continue;
        end
        
        H = H(1:NoSatsUsed, :);
        y = y(1:NoSatsUsed, :);
        
        % W 계산
        W = zeros(NoSatsUsed,NoSatsUsed);
        
        % 고도각 가중치 
        W_el = zeros(NoSatsUsed,NoSatsUsed);
        W_el = diag(sind(matrix_el(:,1)));
        
        % SNR 가중치
        W_snr = zeros(NoSatsUsed,NoSatsUsed);
        
        % 가충치 합산
        W = W_el;

        % 초기화
        matrix_el = 0;

        % xhat 계산
        xhat = pinv(H'*W*H)*H'*W*y;
        % xhat = pinv(H) * y;
        
        % 업데이트
        x = x + xhat;
    
        if norm(xhat) < EpsStop
            % fprintf("[%d] Err : %3.1fm\n", gs, norm(x(1:3) - TruePos(:)));

            nEst = nEst + 1;
            estm(nEst,1) = gs;
            estm(nEst,2:5) = x;
            estm(nEst, 6) = NoSats;
            estm(nEst, 7) = NoSatsUsed;
            break
        end
    end
end
        
%% RMSE Calc

TTs = estm(:, 1);
XYZ = estm(:, 2:4);
NEV = xyz2topo2(XYZ, TruePos);

[rmse, horErr, ~, dim3Err] = nev2rmse(NEV);

%% Figure

close all; clc;
PlotRMSE(TTs, NEV, estm(:,6), estm(:,7));

%% Console disp

fprintf('%-15s : %6.3f [m]\n', 'Horizontal RMSE', rmse(1));
fprintf('%-15s : %6.3f [m]\n', 'Vertical RMSE', rmse(2));
fprintf('%-15s : %6.3f [m]\n\n', '3D RMSE', rmse(3));

fprintf('%-15s : %6.3f [m]\n', 'Max 2D Error', max(horErr));
fprintf('%-15s : %6.3f [m]\n', 'Max 3D Error', max(dim3Err));