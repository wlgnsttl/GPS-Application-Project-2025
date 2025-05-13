function [] = PlotQM(arrQM, prn, obsType)

temp = arrQM(arrQM(:,2) == prn & arrQM(:,3) == obsType, [1,4:6]);

ObsTimeUTC= temp(:,1);
PRObs = temp(:,2);
CPObs = temp(:,3);
DOPObs = temp(:,4);

ObsTimeUTC = gs2gdhms(ObsTimeUTC);
ObsTimeUTC = dms2deg(ObsTimeUTC(:,2),ObsTimeUTC(:,3),ObsTimeUTC(:,4));

figure;

subplot(3,1,1);
hold on; grid on;
xlim([0,24]);
ylabel('Code-Pseudo Range [m]');
plot(ObsTimeUTC,PRObs);
hold off;

subplot(3,1,2);
hold on; grid on;
xlim([0,24]);
ylabel('Carrier-Phase [cycle]');
plot(ObsTimeUTC,CPObs);
hold off;

subplot(3,1,3);
hold on; grid on;
xlim([0,24]);
xlabel('hours');
ylabel('Doppler [Hz]');
plot(ObsTimeUTC,DOPObs);
hold off;

sgtitle(sprintf('PRN[%s], ObsType[%s]', num2str(prn), num2str(obsType)));