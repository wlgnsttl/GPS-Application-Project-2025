function [] = CompEPH(arrQM,sp3, eph, prn)

prn_arrQM = arrQM(arrQM(:,2)==prn, :);
TTs_arrQM = unique(prn_arrQM(:,1));

prn_sp3 = sp3(sp3(:,2)==prn, :);
TTs_sp3 = unique(prn_sp3(:,1));

TTs_diff_arrQM = diff(TTs_arrQM);
interval = mode(TTs_diff_arrQM);

diff_idx = find(TTs_diff_arrQM > interval);

TTs_bounds = zeros(length(diff_idx)+1, 2);
TTs_bounds(2:end,1) = diff_idx+1;
TTs_bounds(1:end-1,2) = diff_idx;
TTs_bounds(1,1) = 1; TTs_bounds(end,2) = length(TTs_arrQM);
TTs_bounds = TTs_arrQM(TTs_bounds);
TTs_bounds = reshape(TTs_bounds, [length(diff_idx)+1, 2]);



TTs = [];
for i = 1:size(TTs_bounds, 1)
    TT_begin = TTs_bounds(i,1);
    TT_end = TTs_bounds(i,2);

    TT_pick = TTs_sp3(TTs_sp3>=TT_begin & TTs_sp3<=TT_end);
    
    TTs = [TTs; TT_pick];
end

xyz_brdc = NaN(length(TTs), 4);
xyz_sp3 = NaN(length(TTs), 4);
xyz_err = NaN(length(TTs), 5);

for epoch = 1:length(TTs)
    ieph = PickEPH_multi(eph, prn, TTs(epoch));
    xyz_brdc(epoch, :) = [TTs(epoch),getSatPos_lab(eph, ieph, TTs(epoch))];
    xyz_sp3(epoch, :) = prn_sp3(prn_sp3(:,1) == TTs(epoch), [1, 3:5]);

    xyz_err(epoch, 1:4) = [TTs(epoch), xyz_brdc(epoch, 2:4) - xyz_sp3(epoch, 2:4)];
    xyz_err(epoch, 5) = norm(xyz_err(epoch, 2:4));

    fprintf("%3d     %6d %14.3f %14.3f %14.3f (BRDC)\n", ...
            prn, TTs(epoch), xyz_brdc(epoch, 2:4));
    
    fprintf("               %14.3f %14.3f %14.3f (SP3)\n", ...
            xyz_sp3(epoch, 2:4));
    
    fprintf("               %14.3f %14.3f %14.3f (3D)  %3.1fm\n\n", ...
            xyz_err(epoch, 2:5));
    
end

TimeUTC = gs2gdhms(TTs);
TimeUTC = dms2deg(TimeUTC(:,2),TimeUTC(:,3),TimeUTC(:,4));

%%
figure;

subplot(4,2,1);
hold on
ylabel('X[m]'); xlim([0, 24]); grid on;

plot(TimeUTC, xyz_sp3(:,2), 'xr');
plot(TimeUTC, xyz_brdc(:,2), ':og');
legend('BRDC', 'SP3');

%%
subplot(4,2,3);
hold on
ylabel('Y[m]'); xlim([0, 24]); grid on;

plot(TimeUTC, xyz_sp3(:,3), 'xr');
plot(TimeUTC, xyz_brdc(:,3), ':og');

%%
subplot(4,2,5);
hold on
ylabel('Z[m]'); xlim([0, 24]); grid on;

plot(TimeUTC, xyz_sp3(:,4), 'xr');
plot(TimeUTC, xyz_brdc(:,4), ':og');
xlabel('Hours');

%%
subplot(4,2,2);
hold on
ylabel('\DeltaX[m]'); xlim([0, 24]); grid on;

scatter(TimeUTC, xyz_err(:,2), 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',0.1);

%%
subplot(4,2,4);
hold on
ylabel('\DeltaY[m]'); xlim([0, 24]); grid on;

scatter(TimeUTC, xyz_err(:,3), 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',0.1);

%%
subplot(4,2,6);
hold on
ylabel('\DeltaZ[m]'); xlim([0, 24]); grid on;

scatter(TimeUTC, xyz_err(:,4), 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',0.1);

%%
subplot(4,2,8);
hold on
ylabel('\Delta 3D[m]'); xlim([0, 24]); grid on;

plot(TimeUTC, xyz_err(:,5),':o', 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',0.1, 'Color', 'k');
xlabel('Hours');




