clear; close all;
clc;

addpath(genpath("data\"));
addpath(genpath("functions\"));

load('QM_RTAP1_250425_0659.mat');
load('eph_25115_1.mat');
load('TruePos_RTAP1_250425_0659.mat');

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
estm = zeros(NoEpochs, 9);
nEst = 0;

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
            H(NoSatsUsed,:) = [g', k', 1];
            y(NoSatsUsed) = obs_dopp - com;
        end
        
        if NoSatsUsed < 7
            continue;
        end

        H = H(1:NoSatsUsed, :);
        y = y(1:NoSatsUsed, :);

        xhat = pinv(H) * y;

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

TTs = estm(:, 1);
XYZ = estm(:, 2:4);
VXYZ = estm(:, 5:7);

NEV = xyz2topo3(XYZ, TruePos);
VNEV = xyz2topo3(VXYZ, TrueVel);

[rmse, horErr, verErr, dim3Err] = nev2rmse(NEV);

%% Figure

close all; clc;
PlotRMSE(TTs, NEV, estm(:,9), estm(:,10));
PlotRMSE(TTs, VNEV, estm(:,9), estm(:,10));
