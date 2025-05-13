clear;
clc;
close all

[arrQM, finalPRNs, finalTTs] = ReadQM('QM_GAMG00KOR_R_20240010000_01D_30S_MO');
%%
PlotQM(arrQM, 502, 103);