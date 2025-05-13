clear
clc
close all

%%
PRN_num = 21;
SYS = 100;

sp3 = ReadSP3('IGS0OPSFIN_20250220000_01D_15M_ORB.SP3');
eph = ReadEPH_multi('BRDC00IGS_R_20250220000_01D_MN.rnx');

%%
sp3_pos = sp3(sp3(:,2) == PRN_num+SYS,:);
TTs = unique(sp3_pos(:,1));

rnx_pos = NaN(length(TTs), 3);

for epoch = 1:length(TTs)
    ieph = PickEPH_multi(eph,PRN_num+SYS,TTs(epoch));
    rnx_pos(epoch, :) = getSatPos_lab(eph, ieph, TTs(epoch));
end

err_val = rnx_pos - sp3_pos(:,3:5);

%%
TimeUTC = gs2gdhms(TTs);
TimeUTC = dms2deg(TimeUTC(:,2),TimeUTC(:,3),TimeUTC(:,4));

%%
figure; clf;

subplot(2,3,1);
hold on
title('X'); xlabel('hours'); ylabel('meters'); xlim([0, 25]); grid on;

plot(TimeUTC, sp3_pos(:,3), 'b', LineWidth=2);
plot(TimeUTC, rnx_pos(:,1), 'g--', LineWidth=2);
legend('SP3', 'RINEX');

%%
subplot(2,3,2);
hold on
title('Y'); xlabel('hours'); ylabel('meters'); xlim([0, 25]); grid on;

plot(TimeUTC, sp3_pos(:,4), 'b', LineWidth=2);
plot(TimeUTC, rnx_pos(:,2), 'g--', LineWidth=2);
legend('SP3', 'RINEX');

%%
subplot(2,3,3);
hold on
title('Z'); xlabel('hours'); ylabel('meters'); xlim([0, 25]); grid on;

plot(TimeUTC, sp3_pos(:,5), 'b', LineWidth=2);
plot(TimeUTC, rnx_pos(:,3), 'g--', LineWidth=2);
legend('SP3', 'RINEX');

%%
subplot(2,3,4);
hold on
title('Orbit Error X'); xlabel('hours'); ylabel('meters'); xlim([0, 25]); grid on;

plot(TimeUTC, err_val(:,1), 'b', LineWidth=2);

%%
subplot(2,3,5);
hold on
title('Orbit Error Y'); xlabel('hours'); ylabel('meters'); xlim([0, 25]); grid on;

plot(TimeUTC, err_val(:,2), 'b', LineWidth=2);

%%
subplot(2,3,6);
hold on
title('Orbit Error Z'); xlabel('hours'); ylabel('meters'); xlim([0, 25]); grid on;

plot(TimeUTC, err_val(:,3), 'b', LineWidth=2);


