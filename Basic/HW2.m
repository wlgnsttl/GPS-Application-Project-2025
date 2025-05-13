clear;
clc;
close all;

%%
[arrQM] = ReadQM('QM_GAMG00KOR_R_20240010000_01D_30S_MO');
sp3 = ReadSP3('COD0OPSFIN_20240010000_01D_05M_ORB.SP3');
eph = ReadEPH_multi('BRDM00DLR_S_20240010000_01D_MN.rnx');

%%
prn = 132;
CompEPH(arrQM, sp3, eph, prn);