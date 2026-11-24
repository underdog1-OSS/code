function estimate = wdr_oshea_estimate(x, time, ridgeTime, ridge, mList, cfg, pairs)
% 先完整运行原 WD-RANSAC，仅对最终最高分直线做一次 O'Shea 精化。
% 然后比较原最高分直线与精化直线的相干得分，保留较高者。
% 得分严格提高才采用精化直线，相等或降低时保留原直线。
% 当前对照只设一个 M，所以每次试验只精化一条直线。
% 若传入多个 M，则每个预算只精化各自原始最高分直线一次。
x = x(:); time = time(:); ridgeTime = ridgeTime(:); ridge = ridge(:);

% 1. 直接复用原基线，完成全部候选的生成、评分和最高分选择。
if nargin<7
    original = wdr_estimate(x,time,ridgeTime,ridge,mList,cfg);
else
    original = wdr_estimate(x,time,ridgeTime,ridge,mList,cfg,pairs);
end

% 2. 只精化原最高分直线；不再对全部抽样候选做精化。
[aFine,bFine] = wdr_oshea_refine( ...
    x,time,original.a,original.b,cfg.osheaWindowLength);

% 3. 用全部原始样本重新计算精化直线的相干得分，与原最高分比较。
% 原直线得分直接沿用原基线已经算出的值，评分定义完全相同。
finePhase = (time.^2/2)*aFine+time*bFine;
fineScores = abs(x.'*exp(-1j*finePhase));
useFine = fineScores>original.score;
estimate = original;
estimate.a(useFine) = aFine(useFine);
estimate.b(useFine) = bFine(useFine);
estimate.score(useFine) = fineScores(useFine);
estimate.omega = ridgeTime*estimate.a+ones(size(ridgeTime))*estimate.b;

% 保留原粗候选列表及其最高分来源；最终输出参数可能已精化。
% bestIteration 和 bestPairs 始终指向原 WD-RANSAC 胜出的候选。
estimate.coarseCandidateA = original.candidateA;
estimate.coarseCandidateB = original.candidateB;
estimate.coarseCandidateScores = original.candidateScores;
estimate.originalBestA = original.a;
estimate.originalBestB = original.b;
estimate.originalBestScore = original.score;
estimate.refinedBestA = aFine;
estimate.refinedBestB = bFine;
estimate.refinedBestScore = fineScores;
estimate.usedRefinement = useFine;
estimate.refinementCount = numel(original.a);
estimate.osheaWindowLength = cfg.osheaWindowLength;
end
