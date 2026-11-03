function estimate = wdr_estimate(x, time, ridgeTime, ridge, mList, cfg, pairs)% 原始信号x，对应时间time，中间256时刻ridgeTime，每个时刻的峰值ridge
% 估计模块：两个峰值点拟合一条直线，用全部 512 点的相干补偿幅度评分。              要比较的选点轮数mList，配置参数cfg，pairs提前指定要选哪些点对
% 输出评分最高的候选直线；对应原文 Step 3--6 与式 (3)。
% 不进行一致点计数、最小二乘重估或局部优化，这些不属于本篇主算法。
%WDR_ESTIMATE Steps 3-6: two-point lines, coherent score, maximum-score line.
% Optional pairs (2-by-max(M)) allows deterministic numerical verification.
% No fit against the truth and no local optimization are performed.
x = x(:); time = time(:); ridgeTime = ridgeTime(:); ridge = ridge(:);
mList = mList(:).';
validateattributes(mList, {'numeric'}, {'integer','positive','finite','nonempty'});% 检查mList是否符合要求
assert(numel(x)==numel(time) && numel(ridgeTime)==numel(ridge), ...                % 检查两组数据的长度是否对应
    'wdr:Dimensions', 'Signal/time or ridge/time lengths do not match.');
maxM = max(mList);                                                                 % 选出最大的选点轮数
if nargin < 7                                                                      % nargin表示传入了几个参数，如果不传pairs（只传6个）
    pairs = wdr_sample_pairs(numel(ridge),maxM,cfg);                               % 从 256 个峰值位置中随机生成 3000 对编号，生成2*3000矩阵
else                                                                               % 检查pairs合规（必须是2*maxM，编号是整数，fix能去掉小数，
    assert(isequal(size(pairs),[2 maxM]) && all(pairs(:)==fix(pairs(:))) && ...    % 编号在1到256之间 
        all(pairs(:)>=1) && all(pairs(:)<=numel(ridge)) && ...                     % 每一列的两个编号不同
        all(pairs(1,:)~=pairs(2,:)), 'wdr:Pairs', 'Invalid time-index pairs.');
end
i1 = pairs(1,:); i2 = pairs(2,:);                                                  % 取出pairs的两行，i1保存每对点的第一个编号，i2保存每对点的第二个编号
a = ((ridge(i2)-ridge(i1))./(ridgeTime(i2)-ridgeTime(i1))).';                      % 计算直线斜率
b = (ridge(i1)-a(:).*ridgeTime(i1)).';                                             % 利用斜率算出截距
scores = zeros(1,maxM);                                                            % 创建一行maxM个零，用来存放每条候选直线的评分
% Score all 512 original complex samples. .'' below is a nonconjugate
% transpose: using x'' would conjugate x and reverse the compensation sign.
for first = 1:cfg.scoreBlockSize:maxM                                              % first是每轮评分的起始直线编号
    block = first:min(first+cfg.scoreBlockSize-1,maxM);                            % block是这一轮要评分的直线编号    
    phase = (time.^2/2)*a(block) + time*b(block);                                  % 计算这轮要评分的直线对应的相位
    scores(block) = abs(x.' * exp(-1j*phase));                                     % 用每条候选直线对应的相位补偿原始信号，把补偿后的 512 个值相加，再取相加结果的幅值作为评分 
end
best = zeros(1,numel(mList));                                                      % mList中存放要比较的直线数
estimate.score = zeros(1,numel(mList));                                            % estimate.score是比较mList个直线中评分最高的分数
for k = 1:numel(mList)                                                             % best是最高评分对应的候选直线编号
    [estimate.score(k), best(k)] = max(scores(1:mList(k)));
end
estimate.a = a(best);                                                              % 保存每次比较的最高评分对应的直线斜率和截距
estimate.b = b(best);
estimate.omega = ridgeTime*estimate.a + ones(size(ridgeTime))*estimate.b;          % 计算选中直线在各个时刻的估计频率
estimate.mList = mList;                                                            % 记录每轮比较了多少条直线
estimate.bestIteration = best;                                                     % 记录每轮比较出的最好的直线编号
estimate.bestPairs = pairs(:,best);                                                % 取出这些直线使用的点对
estimate.candidateA = a;                                                           % 保存所有直线的斜率、截距、分数
estimate.candidateB = b;
estimate.candidateScores = scores;
end

