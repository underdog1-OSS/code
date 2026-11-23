function [ridge, detail] = wdr_local_window_ridge(wd, omega, cfg, mode, initialOmega)
% 三组都返回每列一个频率格点；只改变选峰，不改变 WD 和后端。
% mode 为 A/B/C。C 仅在首列用真值最近格点，之后与 B 完全一致。
omega = omega(:);
assert(size(wd,1)==numel(omega) && all(diff(omega)>0), 'WD 维度或频率轴错误。');
allowed = find(abs(omega)<=cfg.omegaLimit+1e-10);
nTime = size(wd,2);
detail.index = zeros(nTime,1);
detail.lower = nan(nTime,1);
detail.upper = nan(nTime,1);
detail.peakValue = zeros(nTime,1);
if strcmp(mode,'A')
    % 全频峰规则直接调用已有模块。
    [ridge,detail.peakValue] = wdr_ridge(wd(allowed,:),omega(allowed));
    [~,index] = max(wd(allowed,:),[],1);
    detail.index = allowed(index(:));
    return;
end
assert(any(strcmp(mode,{'B','C'})), '方法只能为 A、B 或 C。');
if strcmp(mode,'B')
    [~,localIndex] = max(wd(allowed,1));
    detail.index(1) = allowed(localIndex);
else
    assert(nargin==5 && isfinite(initialOmega), 'C 需要首列的真实频率。');
    [~,localIndex] = min(abs(omega(allowed)-initialOmega));
    detail.index(1) = allowed(localIndex);
end
ridge = zeros(nTime,1);
ridge(1) = omega(detail.index(1));
detail.peakValue(1) = wd(detail.index(1),1);
for k = 2:nTime
    % 唯一的状态就是前一列已选中的频率；不再取全频峰作为中心。
    bounds = sort(cfg.windowFactors*ridge(k-1));
    detail.lower(k) = bounds(1);
    detail.upper(k) = bounds(2);
    % 负的噪声初峰也保留；排序上下界不增加任何纠错机制。
    candidates = allowed(omega(allowed)>=bounds(1)-1e-10 & ...
        omega(allowed)<=bounds(2)+1e-10);
    % 前一列频率本身仍在窗口内，因此不需要回退、重启或补点。
    assert(~isempty(candidates), '局部窗口没有格点，不能偷偷回退为全频取峰。');
    [detail.peakValue(k),localIndex] = max(wd(candidates,k));
    detail.index(k) = candidates(localIndex);
    ridge(k) = omega(detail.index(k));
end
end
