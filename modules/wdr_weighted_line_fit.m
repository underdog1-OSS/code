function [a,b,diagnostic] = wdr_weighted_line_fit(time,ridge,weights)
% 输入峰点的频率方向 WLS：min sum(w_i*(ridge_i-a*time_i-b)^2)。
% 零权重点不参与；时间中心化、缩放后用 sqrt(w) 加权的线性求解。
% 不使用法方程，不把权重误平方，不计算额外的相干得分 J。
time = time(:); ridge = ridge(:); weights = weights(:);
assert(numel(time)==numel(ridge) && numel(time)==numel(weights), ...
    'wdr:WeightedFitDimensions','时间、频率和权重长度必须一致。');
validateattributes(time,{'numeric'},{'real','finite','nonempty'});
validateattributes(ridge,{'numeric'},{'real','finite','nonempty'});
validateattributes(weights,{'numeric'},{'real','finite','nonnegative','nonempty'});
active = weights>0;
diagnostic.applied = false;
diagnostic.reason = '正权重点不足两个不同的时刻。';
diagnostic.positivePointCount = sum(active);
a = NaN; b = NaN;
if sum(active)<2 || numel(unique(time(active)))<2, return; end
t = time(active); y = ridge(active);
w = weights(active)/max(weights(active)); % 公共比例不改变 WLS 解。
offset = sum(w.*t)/sum(w);
centeredTime = t-offset;
timeScale = sqrt(sum(w.*centeredTime.^2)/sum(w));
if ~isfinite(timeScale) || timeScale<=0, return; end
rootWeight = sqrt(w);
design = [centeredTime/timeScale ones(size(t))].*rootWeight;
if rank(design)<2, return; end
coefficient = design\(y.*rootWeight);
a = coefficient(1)/timeScale;
b = coefficient(2)-a*offset;
if ~isfinite(a) || ~isfinite(b)
    a = NaN; b = NaN;
    diagnostic.reason = '加权拟合结果不是有限数。';
    return;
end
diagnostic.applied = true;
diagnostic.reason = '使用输入峰点的累计权重完成一次 WLS。';
end
