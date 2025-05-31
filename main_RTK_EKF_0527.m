clc;
clear all;
%% Read Data File
%%: Base (DSAC_17025)
arrQM_base = load('QDSAC_17025');
%---:베이스 참값
TruePos_base    = [-3041235.578 4053941.677 3859881.013];
TrueLLH_base    = xyz2gd(TruePos_base);
LatLon_base_rad = [deg2rad(TrueLLH_base(1)), deg2rad(TrueLLH_base(2))]; 

%%: Rover (DSBC_17025)
arrQM_rov = load('QDSBC_17025');
%---:로버 참값
TruePos_rov    = [-3041241.741 4053944.143 3859873.640];
TrueLLH_rov    = xyz2gd(TruePos_rov);
LatLon_rov_rad = [deg2rad(TrueLLH_rov(1)), deg2rad(TrueLLH_rov(2))];

%%: eph 
eph = ReadEPH_multi('BRDC17025.rnx');
%% 상수 정의
CCC     = 299792458;
obsType = 120; % L1: 120
lambda = CCC/1575420000; % CCC/L1 밴드의 주파수

P = eye(8);
P_code = eye(3)*100;
P_CP = eye(5);
P(1:3,1:3) = P_code;
P(4:8,4:8)= P_CP;

Q = blkdiag(eye(3)*1e-6,eye(5)*1e-12);
% Q = eye(8)*1e-6;
R = eye(10);
%% 저장 변수
% DOP: 
% gs(1) -- # of visible satellite(2) -- PDOP(3) -- HDOP(4) -- VDOP(5)
DOP = [];
nEst = 0;
%% QM 선별 (120 = code PR, 111 = Carrier Phase)
QM_base = arrQM_base(arrQM_base(:,3) == 120 | arrQM_base(:,3) == 111, :);
QM_rov  = arrQM_rov(arrQM_rov(:,3) == 120 | arrQM_rov(:,3) == 111, :);
QM_base = sortrows(QM_base, [1 2]);
QM_rov  = sortrows(QM_rov, [1 2]);

TTs = unique(QM_rov(:,1));

NoSats = 6;
% 모호정수 초기값 10
x = [TruePos_rov, 10, 10, 10, 10, 10];
AppPos = x;
LatLon_AppPos     = xyz2gd(AppPos(1:3));
LatLon_AppPos_rad = [deg2rad(LatLon_AppPos(1)), deg2rad(LatLon_AppPos(2))];

SatVec = [10 12 29 14 31 32];

tts = unique(QM_rov(:,1));
new_QM_rover = [];
for k = 1 : length(tts)
    idx_tts = QM_rov(:,1) == tts(k);
    cQM = QM_rov(idx_tts,:);
    for n = 1 : length(SatVec)
        idx_prn = cQM(:,2) == SatVec(n);
        new_QM_rover = [new_QM_rover; cQM(idx_prn,:)];
    end
end

tts = unique(QM_base(:,1));
new_QM_base = [];
for k = 1 : length(tts)
    idx_tts = QM_base(:,1) == tts(k);
    cQM = QM_base(idx_tts,:);
    for n = 1 : length(SatVec)
        idx_prn = cQM(:,2) == SatVec(n);
        new_QM_base = [new_QM_base; cQM(idx_prn,:)];
    end
end

y_p = [];y_c = [];
H_p = [];H_c = [];

N_total = [];

for kE = 1:length(TTs)-10 %pivot위성 32 유지하는 시간대로만.
% for kE = 1:length(TTs)
    AppPos = x(1:3);
    QM1e_rov = new_QM_rover(new_QM_rover(:,1) == TTs(kE), : );
    QM1e_base = new_QM_base(new_QM_base(:,1) == TTs(kE), : );
    gs     = QM1e_base(1,1);
    %: st base elevation(1) -- st rov elevation(2)
    el1e = zeros(NoSats,2);
    %: st base com(1) -- st rov com(2)
    com = zeros(NoSats,2);
    %: st base xyz(1,2,3) --- st rov xyz(4,5,6)
    Rec2Sat = zeros(NoSats,2);

    %: gs(1) -- PRN(2) -- Band/Channel(3) -- Pseudorange(관측)(4) -- Com_Pseudorange(계산된 거리 by 코드의사거리)(5)
    %: carrierphase * lambda(관측)(6) -- Com_carrierphase(계산된 거리 by 반송파위상)(7)
    %: 수신기 to 위성 (x,y,z)(8,9,10) -- Elevation(11) -- C/No(12)
    % C/No는 RINEX 2.0이라 없다.
    Pivot_r = zeros(1,12);
    Pivot_b = zeros(1,12);
    Other_r = zeros(NoSats-1,12);
    Other_b = zeros(NoSats-1,12);
    
    % 코드의사거리와 반송파위상를 분리해서(코드: 120, 반송파:111) 해야함 
    QM1e_base_p = QM1e_base(QM1e_base(:,3) == 120,:);
    QM1e_base_c = QM1e_base(QM1e_base(:,3) == 111,:);

    QM1e_rov_p  = QM1e_rov(QM1e_rov(:,3) == 120,:);
    QM1e_rov_c  = QM1e_rov(QM1e_rov(:,3) == 111,:);

    for kS = 1:NoSats
        % base
        prn_base  = 100 + QM1e_base_p(kS,2); % RINEX 2.0이므로 100을 더해야 함
        obs_base  = QM1e_base_p(kS,4);
        STT_base  = obs_base/CCC;
        tc_base   = gs - STT_base;
        ieph_base = PickEPH_multi(eph, prn_base, tc_base);
        tmp_base  = getSatPos(eph, ieph_base, tc_base);
        tmp_base  = RotSatPos(tmp_base, STT_base);
        com(kS,1) = norm(tmp_base - TruePos_base);
        Rec2Sat(kS,1:3) = tmp_base - TruePos_base;

        % rov
        prn_rov   = 100 + QM1e_rov_p(kS,2); % RINEX 2.0이므로 100을 더해야 함
        obs_rov   = QM1e_rov_p(kS,4);
        STT_rov   = obs_rov/CCC;
        tc_rov    = gs - STT_rov;
        ieph_rov  = PickEPH_multi(eph, prn_rov, tc_rov);
        tmp_rov   = getSatPos(eph, ieph_rov, tc_rov);
        tmp_rov   = RotSatPos(tmp_rov, STT_rov);
        com(kS,2) = norm(tmp_rov - AppPos(1:3));
        Rec2Sat(kS,4:6) = tmp_rov - AppPos(1:3);

        % el
        el1e(kS,1) = xyz2el(tmp_base - TruePos_base,LatLon_base_rad(1),LatLon_base_rad(2));
        el1e(kS,2) = xyz2el(tmp_rov - AppPos(1:3),LatLon_AppPos_rad(1),LatLon_AppPos_rad(2));
    end
    
    idx_pivot = el1e(:,2) == max(el1e(:,2));
    idx_other = ~idx_pivot;

    if kE == 1
        pivot_st = QM1e_rov_p(idx_pivot,2);
        fprintf('피벗위성 %d\n',pivot_st);
    else
        if pivot_st ~= QM1e_rov_p(idx_pivot,2)
            fprintf('피벗위성 달라짐\n');
            fprintf('%d 이전 피벗: %d, 현재 피벗: %d\n',gs,pivot_st,QM1e_rov_p(idx_pivot,2));
            pivot_st = QM1e_rov_p(idx_pivot,2);
        end
    end
    
    Pivot_b(1,1:4)  = QM1e_base_p(idx_pivot,1:4);
    Pivot_b(1,5)    = com(idx_pivot,1);
    Pivot_b(1,6)    = QM1e_base_c(idx_pivot,4) * lambda;
    % Pivot_b(1,7)    = com(idx_pivot,1) + N(idx_pivot,1) * lambda; % -> 의문 발생, 이 모호정수는 rov의 모호정수이지 base꺼가 아님
    Pivot_b(1,8:10) = Rec2Sat(idx_pivot,1:3);                     % -> 여기에서 사용하는 모호정수는 이중차분된 모호정수이기 때문에 상관X
    Pivot_b(1,11)   = el1e(idx_pivot,1);                          % -> 그냥 여기에서 필요가 없음 밑에서 계산해야 함 
    % Pivot_b(1,12)   = QM1e_base(idx_pivot,7);

    Other_b(:,1:4)  = QM1e_base_p(idx_other,1:4);
    Other_b(:,5)    = com(idx_other,1);
    Other_b(:,6)    = QM1e_base_c(idx_other,4) * lambda;
    % Other_b(:,7)    = com(idx_other,1) + N(idx_other,:) * lambda;
    Other_b(:,8:10) = Rec2Sat(idx_other,1:3);
    Other_b(:,11)   = el1e(idx_other,1);
    % Other_b(:,12)   = QM1e_base(idx_other,7);

    Pivot_r(1,1:4)  = QM1e_rov_p(idx_pivot,1:4);
    Pivot_r(1,5)    = com(idx_pivot,2);
    Pivot_r(1,6)    = QM1e_rov_c(idx_pivot,4) * lambda;
    % Pivot_r(1,7)    = com(idx_pivot,2) + N(idx_pivot,1)*lambda;
    Pivot_r(1,8:10) = Rec2Sat(idx_pivot,4:6);
    Pivot_r(1,11)   = el1e(idx_pivot,2);
    % Pivot_r(1,12)   = QM1e_rov(idx_pivot,7);

    Other_r(:,1:4)  = QM1e_rov_p(idx_other,1:4);
    Other_r(:,5)    = com(idx_other,2);
    Other_r(:,6)    = QM1e_rov_c(idx_other,4) * lambda;
    % Other_r(:,7)    = com(idx_other,2) + N(idx_other,:)*lambda;
    Other_r(:,8:10) = Rec2Sat(idx_other,4:6);
    Other_r(:,11)   = el1e(idx_other,2);
    % Other_r(:,12)   = QM1e_rov(idx_other,7);

    %% obs dd
    % 수신기간 단일차분 먼저 수행 후 위성간 단일차분
    % 코드의사거리, 반송파 위상 각각 
    dd_obs_p = (Pivot_b(1,4) - Pivot_r(1,4)) - (Other_b(:,4) - Other_r(:,4));
    dd_obs_c = (Pivot_b(1,6) - Pivot_r(1,6)) - (Other_b(:,6) - Other_r(:,6));

    %% com dd
    % 코드의사거리, 반송파 위상 각각 
    dd_com_p = (Pivot_b(1,5) - Pivot_r(1,5)) - (Other_b(:,5) - Other_r(:,5));
    N = round((dd_obs_c - dd_obs_p)/lambda);
    N_total = [N_total;N];
    dd_com_c = (Pivot_b(1,5) - Pivot_r(1,5)) - (Other_b(:,5) - Other_r(:,5)) + lambda*N;
    y_temp_p = dd_obs_p - dd_com_p;
    y_temp_c = dd_obs_c - dd_com_c;

    y_p = y_temp_p;
    y_c = y_temp_c;

    %% H
    %열: rov의 xyz(1,2,3) -- 모호정수(4,5,6,7,8,9,10,11,12)
    %행: code(n-1) --- carrier(n-1)
    % 코드의사거리, 반송파위상 두 경우 1,2,3열은 동일, but 모호정수는 차이존재 
    H_Pivot = Pivot_r(1,8:10)  / Pivot_r(1,5);
    H_Other = Other_r(:,8:10) ./ Other_r(:,5);
    
    H_temp  = H_Pivot - H_Other;
    H_temp_p  = [H_temp zeros(NoSats-1,NoSats-1)];
    H_temp_c  = [H_temp lambda * eye(NoSats-1)];

    H_p = H_temp_p;
    H_c = H_temp_c;
    

    El_Pivot_weight =  sind(Pivot_r(:,11));
    El_Othter_weight = sind(Other_r(:,11));

    % R 행렬 초기화
    % R_Other_temp = diag(1 ./ (El_Othter_weight .^ 2)); 
    % R_Pivot_temp = diag(1 ./ (El_Pivot_weight .^ 2));
    R_Other_temp = diag(1 ./ (El_Othter_weight.^2)); 
    R_Pivot_temp = diag(1 ./ (El_Pivot_weight.^2));
    
    R(1:5,1:5) = R_Other_temp*100;
    R(1:5,6:10) = R_Pivot_temp*10;
    R(6:10,6:10) = R_Other_temp;
    R(6:10,1:5) = R_Pivot_temp*10;

    H = [H_p;H_c];
    y = [y_p;y_c];

    K = (P*H')/(H*P*H' + R);
    x = x + (K*y)';
    % P = P + Q;
    P = (eye(8) - K * H) * P + Q;
    
    nEst           = nEst + 1;
    estm(nEst,1)   = gs;
    estm(nEst,2:4) = x(1:3);

    nev_err(nEst,1) = gs;
    nev_err(nEst,2:4) = xyz2topo(TruePos_rov - x(1:3), TrueLLH_rov(1), TrueLLH_rov(2));

end


%% 오차 및 RMSE
% 전체 수평 및 수직 오차 계산
% rmse_n = sqrt(mean(nev_err(:,3).^2));
% horizontal_error = sqrt(nev_err(:,2).^2 + nev_err(:,3).^2);
% vertical_error = abs(nev_err(:,4));
% threed_error = sqrt(nev_err(:,2).^2 + nev_err(:,3).^2 + nev_err(:,4).^2);
% 
% % 2D, 3D RMSE 계산
% rmse_horizontal = sqrt(mean(horizontal_error.^2));
% rmse_vertical = sqrt(mean(vertical_error.^2));
% fprintf('RMSE (Horizontal): %.4f m\n', rmse_horizontal);
% fprintf('RMSE (Vertical): %.4f m\n', rmse_vertical);
% 
% rmse_3d = sqrt(mean(threed_error.^2));
% fprintf('RMSE (3D): %.4f m\n', rmse_3d);
% time_h = size(nev_err);

% scatter(E-N)
% 10cm 원으로 그리기 -> 수cm로 설명하기 편함 
figure('Color','white');
r = 0.5;
theta = linspace(0, 2*pi, 100);  % 각도
x_circle = r * cos(theta);
y_circle = r * sin(theta);
plot(x_circle, y_circle, 'b--', 'LineWidth', 1); % 점선 원
hold on;
plot(nev_err(:,3),nev_err(:,2),'b.','MarkerSize',12); axis([-0.3 0.3 -0.3 0.3]); axis square; hold on;
plot(nev_err(1,3),nev_err(1,2),'r^','MarkerSize',6,'MarkerFaceColor','r');
plot(nev_err(10,3),nev_err(10,2),'m^','MarkerSize',6,'MarkerFaceColor','m');
plot(nev_err(300,3),nev_err(300,2),'g^','MarkerSize',6,'MarkerFaceColor','c');
plot(nev_err(3590,3),nev_err(3590,2),'g^','MarkerSize',6,'MarkerFaceColor','g');
legend('','','1 sec','30 sec','300 sec','3590 sec','autoupdate','off');
xline(0,'k--'); yline(0,'r--'); set(gca,'fontname','times new roman','fontsize',13); ylabel('\DeltaN'); xlabel('\DeltaE');

fprintf('1초 : n방향 %.2fcm, e방향 %.2fcm, v방향 %.2fcm \n',100*nev_err(1,2),100*nev_err(1,3),100*nev_err(1,4));
fprintf('30초 : n방향 %.2fcm, e방향 %.2fcm, v방향 %.2fcm \n',100*nev_err(30,2),100*nev_err(30,3),100*nev_err(30,4));
fprintf('139초 : n방향 %.2fcm, e방향 %.2fcm, v방향 %.2fcm \n',100*nev_err(139,2),100*nev_err(139,3),100*nev_err(139,4));
fprintf('3590초 : n방향 %.2fcm, e방향 %.2fcm, v방향 %.2fcm \n',100*nev_err(3590,2),100*nev_err(3590,3),100*nev_err(3590,4));