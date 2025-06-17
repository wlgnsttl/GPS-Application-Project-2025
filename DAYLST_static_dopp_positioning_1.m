clear; close all;
clc;

addpath(genpath("data/"));
addpath(genpath("functions/"));

% 데이터 로드
load('QM_GAMGKOR.mat');
load('eph_24.mat');
load('TruePos_GAMGKOR.mat');

%% 상수, 변수 정의

CCC = 299792458;
L1_lamda = 0.19029;

sys = [100 400];
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
estm = zeros(1, 8);
NoEpochs = length(FinalTTs);
MaxSnr = max(QM(:,7));

for kIter = 1:MaxIter
    HTH = zeros(7,7);
    Hty = zeros(7,1);
    
    vec_rec_p = x(1:3);
    vec_rec_v = x(4:6);

    for kE = 1:NoEpochs
        idx   = QM(:,1)==FinalTTs(kE);
        QM1e  = QM(idx ,:);
        NoSats= size(QM1e,1);
        gs    = QM1e(1,1);
           
        for kS = 1:NoSats
            prn = QM1e(kS,2);
            obs_dopp = QM1e(kS,6) * -L1_lamda;
    
            ieph  = PickEPH_multi(eph,prn,gs);
    
            if ieph == 0
                continue;
            end
            
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
            
            [~, el] = xyz2azel(vec_rho_p, Truellh(1), Truellh(2));

            % 고도각 가중치 
            W_el = 1;
            W_el = WeightEl(rad2deg(el));
    
            % SNR 가중치
            W_snr = 1;
            W_snr = WeightSNR(QM1e(kS,7),MaxSnr);
            
            % 가충치 합산
            W = W_el .* W_snr;

            H = [g', k', 1];
            y = obs_dopp - com;

            HTH = HTH + H'*W*H;
            Hty = Hty + H'*W*y;
        end
    end

    xhat = HTH \ Hty;
    x = x + xhat;

    if norm(xhat) < EpsStop
        estm(1,1)   = 0;
        estm(1,2:8) = x;
        break
    end
end

%% RMSE Calc

XYZ = estm(2:4);
VXYZ = estm(5:7);
dclk = estm(8);

NEV = xyz2topo2(XYZ, TruePos);
VNEV = xyz2topo2(VXYZ + TruePos, TruePos);
