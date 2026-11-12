function [estimate, quality, schedule, timing] = wdr_prosac_estimate( ...
    x, time, ridgeTime, ridge, wd, omega, mList, cfg, growthBudget)
%WDR_PROSAC_ESTIMATE 保留 WD-RANSAC 评分，改用质量排序及渐进抽样。
% 不使用 PROSAC 的停止规则；每种 M 仍取前 M 条候选直线的最高分。
validateattributes(mList, {'numeric'}, {'vector','integer','positive','finite','nonempty'});
timer = tic;
[~,peakIndex] = max(wd,[],1);
assert(numel(ridge)==size(wd,2) && numel(omega)==size(wd,1), ...
    'wdr:ProsacDimensions', 'WD, ridge and frequency grid dimensions must match.');
assert(isequal(ridge(:),reshape(omega(peakIndex),[],1)), ...
    'wdr:ProsacRidge', 'Ridge must contain the maximum WD frequency at each time.');
options = wdr_quality_config(wd,omega);
quality = wdr_peak_quality(wd,ridgeTime,omega,peakIndex,options);
timing.qualitySeconds = toc(timer);
timer = tic;
[pairs,schedule] = wdr_prosac_pairs(quality.order,max(mList),growthBudget);
timing.samplingSeconds = toc(timer);
timer = tic;
estimate = wdr_estimate(x,time,ridgeTime,ridge,mList,cfg,pairs);
timing.scoringSeconds = toc(timer);
schedule.pairs = pairs;
end
