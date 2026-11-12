function verify_wdr_prosac
%VERIFY_WDR_PROSAC 检查质量分数、渐进抽样和基线兼容性。
% 在 MATLAB 编辑器中点击“运行”；不生成结果文件。
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(scriptDir),'modules'));

% 一条局部线性轨迹中放入一个孤立错误峰，检查留一一致性。
nTime = 12; nFreq = 64;
time = (0:nTime-1)'; omega = (0:nFreq-1)'*pi;
wd = repmat(1+0.05*mod((1:nFreq)',5),1,nTime);
peakIndex = (20+(1:nTime))';
peakIndex(6) = 40;
for i = 1:nTime
    wd(peakIndex(i),i) = 10;
end
options = wdr_quality_config(wd,omega);
quality = wdr_peak_quality(wd,time,omega,peakIndex,options);
assert(all(quality.localUsed));
assert(quality.consistency(6)<0.1);
assert(quality.q(6)<median(quality.q));
assert(isequal(sort(quality.order),(1:nTime)'));

% 边界背景不足时回退至全部可用格点；局部邻居不足时显式标记。
smallWd = zeros(12,4); smallTime = (1:4)'; smallOmega = (1:12)';
smallOptions = wdr_quality_config(smallWd,smallOmega);
smallOptions.backgroundRadiusBins = 4;
smallQuality = wdr_peak_quality(smallWd,smallTime,smallOmega,ones(4,1),smallOptions);
assert(all(smallQuality.backgroundFallback));
assert(all(~smallQuality.localUsed));
assert(all(smallQuality.q==0));
assert(isequal(smallQuality.order,(1:4)'));
smallOptions.minBackgroundBins = 10;
raised = false;
try
    wdr_peak_quality(smallWd,smallTime,smallOmega,ones(4,1),smallOptions);
catch err
    raised = strcmp(err.identifier,'wdr:QualityBackground');
end
assert(raised);

% 检查排序编号映射、两点不同、日程单调和新点抽样，以及可复现性。
ranking = [5 4 3 2 1];
rng(19,'twister'); [pairs,schedule] = wdr_prosac_pairs(ranking,12,5);
rng(19,'twister'); [pairs2,schedule2] = wdr_prosac_pairs(ranking,12,5);
assert(isequal(pairs,pairs2) && isequal(schedule.prefixSize,schedule2.prefixSize));
assert(all(pairs(1,:)~=pairs(2,:)) && all(diff(schedule.prefixSize)>=0));
assert(schedule.prefixSize(1)==2 && schedule.prefixSize(end)==5);
rankOf = zeros(5,1); rankOf(ranking) = 1:5;
for k = 1:size(pairs,2)
    n = schedule.prefixSize(k);
    assert(all(rankOf(pairs(:,k))<=n));
    if schedule.forcedNewest(k)
        assert(any(pairs(:,k)==ranking(n)));
    end
end
assert(~schedule.forcedNewest(end));

% 对照程序的高斯基线与原模块使用相同信号、随机序列和相干评分。
cfg = wdr_config('quick');
cfg.trials = 2; cfg.snrDb = [-5 5]; cfg.mList = [3 10]; cfg.showProgress = false;
old = wdr_monte_carlo(cfg);
comparison = wdr_compare_prosac(cfg,max(cfg.mList));
assert(isequal(old.ransacError,comparison.gaussianError));
assert(all(isfinite(comparison.prosacMse(:))));
assert(all(comparison.correctTop10(:)>=0 & comparison.correctTop10(:)<=1));
assert(abs(comparison.prosacMse(1,1)-mean(comparison.prosacError(:,1,1)))<1e-9);
fprintf('PASS: quality, boundary handling, PROSAC schedule and baseline comparison.\n');
end
