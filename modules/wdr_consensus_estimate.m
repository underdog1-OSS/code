function estimate = wdr_consensus_estimate(x,time,ridgeTime,ridge,cfg,prefixPairs)
% 一致集反馈抽样和末尾历史权重 WLS：总共仍只评价 M 个随机候选。
% 最佳模型只由原始全信号的相干得分 J 决定，不按一致点数量选择。
% B 次初始化结束后计分一次；以后仅在最佳 J 刷新时计分。
% 权重不改变抽样或固定的一致集门限；最终仅在搜索结束时的局部集内重估。
% 最终 WLS 额外计算一次全信号 J，仅严格优于搜索最佳模型时接受。
% 随机候选仍为 M 个；有效 WLS 使总 J 评价次数增加一次。
% cfg.consensusTau 必须由实际频率网格计算，不是以 Hz 或编号直接代入。
x = x(:); time = time(:); ridgeTime = ridgeTime(:); ridge = ridge(:);
M = cfg.mList; B = cfg.initialBudget; tau = cfg.consensusTau;
validateattributes(M,{'numeric'},{'scalar','integer','positive'});
validateattributes(B,{'numeric'},{'scalar','integer','positive','<',M});
validateattributes(tau,{'numeric'},{'scalar','real','finite','positive'});
assert(numel(x)==numel(time) && numel(ridge)==numel(ridgeTime), ...
    '信号、时间轴和峰点的维度不一致。');
assert(numel(ridge)>=2 && all(diff(ridgeTime)>0) && all(isfinite(ridge)), ...
    '需要至少两个不同时刻的有限峰点。');
if nargin<6
    prefixPairs = wdr_sample_pairs(numel(ridge),B,cfg);
end
% 这里确实评价 B 个候选，直接复用原模块的分块评分，不是刷新最佳值 B 次。
initial = wdr_estimate(x,time,ridgeTime,ridge,B,cfg,prefixPairs);
globalA = initial.a; globalB = initial.b; globalScore = initial.score;
bestIteration = initial.bestIteration;
residual = abs(ridge-(globalA*ridgeTime+globalB));
members = find(residual<=tau);
initialMembers = members;
[weightState,increment,thresholds] = wdr_consensus_weight_update([],residual,tau);
weightUpdateIterations = B;
weightUpdateBestIterations = initial.bestIteration;
weightThresholdHistory = thresholds(:);
weightMaxResidualHistory = weightState.lastMaxResidual;
weightIncrementHistory = increment;

candidateA = zeros(1,M); candidateB = zeros(1,M); scores = zeros(1,M);
candidatePairs = zeros(2,M);
candidateA(1:B) = initial.candidateA;
candidateB(1:B) = initial.candidateB;
scores(1:B) = initial.candidateScores;
candidatePairs(:,1:B) = prefixPairs;
bestScoreHistory = zeros(1,M);
bestScoreHistory(1:B) = cummax(scores(1:B));
setSizeBeforeSampling = zeros(1,M-B);
fallbackFlags = false(1,M-B);
updateFlags = false(1,M-B);
halfTimeSquared = time.^2/2;

% 后续抽样依赖每次即时更新的集合，因此逐个候选处理，不提前批量抽取。
for iteration = B+1:M
    step = iteration-B;
    nLocal = numel(members);
    setSizeBeforeSampling(step) = nLocal;
    if nLocal>=2
        % 从集合中均匀抽两个不同元素；跨候选可以重复抽到同一对点。
        first = randi(nLocal);
        second = randi(nLocal-1);
        second = second+(second>=first);
        pair = sort(members([first second]));
    else
        % 仅当不足两个有效时刻时回退全局；仍只生成和评价一个候选。
        pair = wdr_sample_pairs(numel(ridge),1,cfg);
        fallbackFlags(step) = true;
    end
    i = pair(1); j = pair(2);
    a = (ridge(j)-ridge(i))/(ridgeTime(j)-ridgeTime(i));
    b = ridge(i)-a*ridgeTime(i);
    phase = halfTimeSquared*a+time*b;
    score = abs(x.'*exp(-1j*phase));
    candidateA(iteration) = a;
    candidateB(iteration) = b;
    scores(iteration) = score;
    candidatePairs(:,iteration) = pair(:);
    if score>globalScore
        globalA = a; globalB = b; globalScore = score;
        bestIteration = iteration;
        updateFlags(step) = true;
        % 必须用新最佳模型重新筛全部原始峰点，允许旧点退出、新点进入。
        residual = abs(ridge-(globalA*ridgeTime+globalB));
        members = find(residual<=tau); % 此门限始终为 tau，不随计分门限改变。
        [weightState,increment,thresholds] = ...
            wdr_consensus_weight_update(weightState,residual,tau);
        weightUpdateIterations(end+1) = iteration; %#ok<AGROW>
        weightUpdateBestIterations(end+1) = iteration; %#ok<AGROW>
        weightThresholdHistory(:,end+1) = thresholds(:); %#ok<AGROW>
        weightMaxResidualHistory(end+1) = weightState.lastMaxResidual; %#ok<AGROW>
        weightIncrementHistory(:,end+1) = increment; %#ok<AGROW>
    end
    bestScoreHistory(iteration) = globalScore;
end

% WD 适配：只用搜索最佳模型对应的最终局部集，沿用其中各点的累计权重。
% 全局历史权重保留；最终局部集之外的点不参与末尾 WLS。
weights = weightState.weights;
refitWeights = zeros(size(weights));
refitWeights(members) = weights(members);
[weightedA,weightedB,fitDiagnostic] = wdr_weighted_line_fit( ...
    ridgeTime(members),ridge(members),weights(members));
weightedOmega = ridgeTime*weightedA+weightedB;
weightedScore = NaN;
additionalScoreCount = 0;
weightedRefitAccepted = false;
% “原 J”是本方法完成 M 次搜索后的最佳 J，不是另跑的基线或初始化 J。
finalA = globalA; finalB = globalB; finalScore = globalScore;
if fitDiagnostic.applied
    weightedPhase = halfTimeSquared*weightedA+time*weightedB;
    weightedScore = abs(x.'*exp(-1j*weightedPhase));
    additionalScoreCount = 1;
    if weightedScore>globalScore
        finalA = weightedA; finalB = weightedB; finalScore = weightedScore;
        weightedRefitAccepted = true;
    end
end
% 平分、加权 J 更低或拟合退化均保留搜索模型；真值不参与此选择。
estimate.a = finalA; estimate.b = finalB; estimate.score = finalScore;
estimate.omega = ridgeTime*finalA+finalB;
estimate.searchBestA = globalA; estimate.searchBestB = globalB;
estimate.searchBestScore = globalScore;
estimate.searchBestOmega = ridgeTime*globalA+globalB;
estimate.outputScoreEvaluated = true;
estimate.candidateScoreCount = M;
estimate.additionalScoreCount = additionalScoreCount;
estimate.totalScoreCount = M+additionalScoreCount;
% 保存未经择优的 WLS 候选，供计时外计算逐次误差，不能用最终输出代替它。
estimate.weightedRefitA = weightedA; estimate.weightedRefitB = weightedB;
estimate.weightedRefitOmega = weightedOmega;
estimate.weightedRefitScore = weightedScore;
estimate.finalWeights = weights; % 全部原始点的历史累计权重。
estimate.refitWeights = refitWeights; % WLS 实际使用的权重；最终局部集外为零。
estimate.refitPositiveWeightCount = fitDiagnostic.positivePointCount;
estimate.weightUpdateCount = weightState.updateCount;
estimate.weightUpdateIterations = weightUpdateIterations;
estimate.weightUpdateBestIterations = weightUpdateBestIterations;
estimate.weightThresholdHistory = weightThresholdHistory;
estimate.weightMaxResidualHistory = weightMaxResidualHistory;
estimate.weightIncrementHistory = weightIncrementHistory;
estimate.weightedRefitApplied = fitDiagnostic.applied;
estimate.weightedRefitAccepted = weightedRefitAccepted;
estimate.weightedRefitDiagnostic = fitDiagnostic;
% 三个阶段用同一套局部有效权重统计 SSE，与本次 WLS 的目标保持一致。
estimate.preRefitWeightedSse = sum(refitWeights.*(ridge-globalA*ridgeTime-globalB).^2);
estimate.postRefitWeightedSse = sum(refitWeights.*(ridge-weightedOmega).^2);
estimate.outputWeightedSse = sum(refitWeights.*(ridge-finalA*ridgeTime-finalB).^2);
% 以下仅统计局部集外仍带历史权重的点，它们已被排除在本次 WLS 之外。
estimate.weightedPointsOutsideFinalConsensusCount = sum(weights>0 & ...
    abs(ridge-globalA*ridgeTime-globalB)>tau);
estimate.mList = M;
estimate.bestIteration = bestIteration; % 搜索最佳候选编号，不指末尾 WLS 的来源。
estimate.bestPairs = candidatePairs(:,bestIteration);
estimate.candidateA = candidateA;
estimate.candidateB = candidateB;
estimate.candidateScores = scores;
estimate.candidatePairs = candidatePairs;
estimate.bestScoreHistory = bestScoreHistory;
estimate.candidateCount = M;
estimate.initialCandidateCount = B;
estimate.initialBestIteration = initial.bestIteration;
estimate.initialA = initial.a; estimate.initialB = initial.b;
estimate.initialScore = initial.score;
estimate.initialConsensusIndex = initialMembers;
estimate.finalConsensusIndex = members;
estimate.setSizeBeforeSampling = setSizeBeforeSampling;
estimate.updateFlags = updateFlags;
estimate.updateCount = sum(updateFlags);
estimate.fallbackFlags = fallbackFlags;
estimate.fallbackCount = sum(fallbackFlags);
estimate.localCandidateCount = M-B-estimate.fallbackCount;
estimate.tau = tau;
end
