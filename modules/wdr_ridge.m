function [omegaHat, peakValue] = wdr_ridge(wd, omegaGrid)               % omegaGrid存每一行对应的实际频率
% 峰值模块：每个时刻取 WD 最大值对应的频率，即论文式 (2)。
%WDR_RIDGE Eq.(2): position of the maximum signed WD value at each time.
% Taking abs(wd) would incorrectly allow negative WD lobes to win.
[peakValue, index] = max(wd, [], 1);            % 对wd的每一列求最大值，peakValue是wd峰值，index是峰值所在的行数
omegaHat = omegaGrid(index);                    % 峰值所在行对应的频率
omegaHat = omegaHat(:);                         
peakValue = peakValue(:);
end

