function va = wdr_va_ridge(wd, omega, lambda, epsilonRelative)
% F^RR、eta=0 的精确全局脊线搜索；输入为频率格点×时间列。
% Algorithm 1 使用时间×频率矩阵并最大化得分，本实现等价地最小化代价。
% 递推中保留每个频率终点的最优部分路径，不采用局部峰跟随。
if nargin<4, epsilonRelative = 1e-12; end
validateattributes(wd,{'numeric'},{'real','finite','nonempty','2d'});
validateattributes(omega,{'numeric'},{'real','finite','vector','numel',size(wd,1)});
validateattributes(lambda,{'numeric'},{'real','finite','scalar','nonnegative'});
validateattributes(epsilonRelative,{'numeric'},{'real','finite','scalar','>',0,'<',1});
omega = omega(:); [nFreq,nTime] = size(wd);
assert(nFreq>=2 && all(diff(omega)>0),'频率轴必须包含至少两个递增格点。');
deltaOmega = mean(diff(omega));
assert(max(abs(diff(omega)-deltaOmega))<1e-9*max(1,deltaOmega), ...
    '本实现的格点转移惩罚要求均匀频率网格。');

% 数据项按 Algorithm 1 的 log(TF/sum(TF)) 定义处理。
% 先统一缩放以避免大数溢出，再设相对下限，保持幅度缩放不改变路径。
positive = max(double(wd),0);
scale = max(positive(:));
va.allZeroInput = scale==0;
if scale==0
    normalizedPositive = ones(size(wd));
else
    normalizedPositive = max(positive/scale,epsilonRelative);
end
logMass = log(sum(normalizedPositive(:)));
nodeCost = -log(normalizedPositive)+logMass;
% 每列减去相同常数，只改善数值尺度，不改变任何完整路径之间的代价差。
columnOffset = min(nodeCost,[],1);
nodeCost = bsxfun(@minus,nodeCost,columnOffset);

states = (1:nFreq).';
transition = lambda*(states-states.').^2;
cumulative = zeros(nFreq,nTime);
if nFreq<=65535
    predecessor = zeros(nFreq,nTime,'uint16');
else
    predecessor = zeros(nFreq,nTime,'uint32');
end
% 全部起始频率都作为有效初始状态，不用真值或首列最大峰限定起点。
cumulative(:,1) = nodeCost(:,1);
for n = 2:nTime
    % 行是上一时刻的频率 j，列是当前频率 k；每个 k 比较全部 j。
    candidates = bsxfun(@plus,cumulative(:,n-1),transition);
    [minimum,parent] = min(candidates,[],1);
    cumulative(:,n) = nodeCost(:,n)+minimum(:);
    predecessor(:,n) = parent(:);
end
[minimum,lastState] = min(cumulative(:,end));
index = zeros(nTime,1);
index(end) = lastState;
for n = nTime:-1:2
    index(n-1) = predecessor(index(n),n);
end
va.index = index;
va.omega = omega(index);
va.cumulativeCost = cumulative;
va.predecessor = predecessor;
va.nodeCost = nodeCost;
va.columnOffset = columnOffset;
va.finalCostByState = cumulative(:,end);
va.shiftedTotalCost = minimum;
va.totalCost = minimum+sum(columnOffset);
va.dataCost = sum(nodeCost(sub2ind(size(wd),index,(1:nTime).')))+sum(columnOffset);
va.transitionCost = lambda*sum(diff(index).^2);
va.lambdaBin = lambda;
va.lambdaOmega = lambda/deltaOmega^2;
va.deltaOmega = deltaOmega;
va.epsilonRelative = epsilonRelative;
va.clippedFraction = mean(wd(:)<=0);
end
