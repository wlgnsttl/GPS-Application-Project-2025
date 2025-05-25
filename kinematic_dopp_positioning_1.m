clear; close all;
clc;

addpath(genpath("data/"));
addpath(genpath("functions/"));

load('QM_RTAP1_250425_0659.mat');
load('eph_25115_1.mat');
load('TruePos_RTAP1_250425_0659.mat');

sp3 = ReadSP3('IGS0OPSULT_20251131800_02D_15M_ORB.SP3');

TruePos = TruePos(:,2:4);
TrueVel = [0 0 0; diff(TruePos)];

%% 상수, 변수 정의        

CCC = 299792458;
L1_lamda = 0.19029;

sys = 100;
obsType = 103;

% Truellh = xyz2gd(TruePos);

%% QM 선별

QM = SelectQM(arrQM, sys, obsType);
FinalTTs = unique(QM(:,1));

%% 추정에 필요한 반복 조건 및 초기값 설정

MaxIter = 5;
EpsStop = 1e-5; 
x = [TruePos(1,:) 0 0 0 0]; 
x = x';

%% 추정과정 시작

NoEpochs = length(FinalTTs);
estm = zeros(NoEpochs, 10);
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
            obs = QM1e(kS,4);
            
            STT = obs/CCC;
            
            tc = gs - STT;
    
            ieph  = PickEPH_multi(eph,prn,gs);
            
            if eph(ieph, 19) > 0
                continue;
            end  

            b = eph(ieph, 4);
            
            % brdc기반 속도추정
            % [vec_sat_p, ~] = getSatPos_lab(eph, ieph, gs);
            % vec_sat_p = vec_sat_p';
            % 
            % vec_sat_v = getSatVel(eph, ieph, gs)';
            
            % 미분정의기반 속도추정
            [vec_sat_p, vec_sat_v] = getSatVel_diff(eph, ieph, tc, STT);

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
            llh = xyz2gd(vec_rec_p');
            [~, el] = xyz2azel(vec_rho_p, llh(1), llh(2));
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
        W_snr = WeightSNR(matrix_snr, MaxSnr);
        
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
            estm(kE,1)   = gs;
            estm(kE,2:8) = x;
            llh_f(nEst,1:3) = xyz2gd(x(1:3)');
            estm(kE, 9) = NoSats;
            estm(kE, 10) = NoSatsUsed;
            break
        end
    end
end
%% velocity calc
estm_v = zeros(NoEpochs,4);
estm_v(1,1:4) = estm(1,1:4);
estm_v(:,1) = estm(:,1);
        
%% RMSE Calc

idx = estm(:, 1) ~= 0;
idx_v = estm_v(:, 1) ~= 0;
estm = estm(idx, :); TruePos = TruePos(idx, :); TrueVel = TrueVel(idx, :);
estm_v = estm_v(idx_v,:); TruePos_v = TruePos(idx_v, : );

estm_v(1,2:4) = TruePos_v(1,1:3);
for i = 1 : 412
    estm_v(i+1,1) = estm(i+1,1);
    estm_v(i+1,2:4) = estm_v(i,2:4) + estm(i,5:7);
end
llh_v = xyz2gd(estm_v(:,2:4));

TTs = estm(:, 1);
XYZ = estm(:, 2:4);
VXYZ = estm(:, 5:7);

TTs_v = estm_v(:, 1);
XYZ_v = estm_v(:, 2:4);

NEV = xyz2topo3(XYZ, TruePos);
VNEV = xyz2topo3(VXYZ, TrueVel);

NEV_v = xyz2topo3(XYZ_v, TruePos_v);

[rmse, horErr, verErr, dim3Err] = nev2rmse(NEV);

[rmse_v, horErr_v, verErr_v, dim3Err_v] = nev2rmse(NEV);
%% Figure

close all; clc;
PlotPosRMSE(TTs, NEV, estm(:,9), estm(:,10));
PlotVelRMSE(TTs, VNEV, estm(:,9), estm(:,10));
PlotPosRMSE(TTs, NEV_v, estm(:,9), estm(:,10));

llh_t = xyz2gd(TruePos);
figure;
geoplot(llh_v(:,1),llh_v(:,2),'b');
figure;
geoplot(llh_t(:,1),llh_t(:,2),'r');
%% Console disp

fprintf('%-15s : %6.3f [m]\n', 'Horizontal RMSE', rmse(1));
fprintf('%-15s : %6.3f [m]\n', 'Vertical RMSE', rmse(2));
fprintf('%-15s : %6.3f [m]\n\n', '3D RMSE', rmse(3));

fprintf('%-15s : %6.3f [m]\n', 'Max 2D Error', max(horErr));
fprintf('%-15s : %6.3f [m]\n', 'Max 3D Error', max(dim3Err));