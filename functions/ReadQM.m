function [arrQM, FinalPRNs, FinalTTs] = ReadQM(filepath)

arrQM = load(filepath);
FinalPRNs = unique(arrQM(:,2));
FinalTTs = unique(arrQM(:,1));