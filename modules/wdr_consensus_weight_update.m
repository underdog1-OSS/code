function [state,increment,thresholds] = wdr_consensus_weight_update(state,residual,tau)
% 多阈值历史支持计分。state=[] 时初始化；本函数不使用信号真值或幅度。
% 首次最大门限为 tau；以后用 (tau+上次最宽内点集的最大残差)/2。
% 每通过一个严格小于门限的判定加一分，未通过时不扣分。
residual = residual(:);
validateattributes(residual,{'numeric'},{'real','finite','nonnegative','nonempty'});
validateattributes(tau,{'numeric'},{'scalar','real','finite','positive'});
if isempty(state)
    state.weights = zeros(size(residual));
    state.updateCount = 0;
    state.lastMaxResidual = NaN;
else
    assert(isequal(size(state.weights),size(residual)), ...
        'wdr:WeightDimensions','累计权重必须对应全部原始峰点。');
end
if state.updateCount==0
    outerThreshold = tau;
else
    outerThreshold = (tau+state.lastMaxResidual)/2;
end
thresholds = outerThreshold*(1:3)/3;
increment = sum(residual<thresholds,2);
inside = residual<thresholds(3);
assert(any(inside),'wdr:EmptyWeightSupport', ...
    '当前模型在最宽加权门限内没有支持点，无法定义下次门限。');
state.weights = state.weights+increment;
state.updateCount = state.updateCount+1;
state.lastMaxResidual = max(residual(inside));
state.lastThresholds = thresholds;
end
