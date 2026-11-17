function [cost, predecessor] = wdr_va_recurrence(nodeCost, lambda)
%WDR_VA_RECURRENCE 每个频率终点的最优部分路径，频率 x 时间。
validateattributes(nodeCost,{'numeric'},{'real','finite','nonempty','2d'});
validateattributes(lambda,{'numeric'},{'real','finite','scalar','nonnegative'});
[nFreq,nTime] = size(nodeCost);
indices = (1:nFreq)';
penalty = lambda*bsxfun(@minus,indices,indices').^2;
assert(all(isfinite(penalty(:))), 'wdr:VaOverflow', 'Transition cost overflow.');
cost = zeros(nFreq,nTime);
predecessor = zeros(nFreq,nTime);
cost(:,1) = nodeCost(:,1);
for n = 2:nTime
    % 行是当前终点 k，列是上一时刻前驱 j；首个并列前驱优先。
    alternatives = bsxfun(@plus,penalty,cost(:,n-1)');
    [best,predecessor(:,n)] = min(alternatives,[],2);
    cost(:,n) = nodeCost(:,n)+best;
end
assert(all(isfinite(cost(:))), 'wdr:VaOverflow', 'Cumulative cost overflow.');
end
