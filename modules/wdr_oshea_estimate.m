function estimate = wdr_oshea_estimate(x, time, ridgeTime, ridge, mList, cfg, pairs)
% 每条候选精化后，比较原直线与精化直线的相干得分，保留较高者。
% 原 wdr_estimate.m 保持独立，作为不精化的基线。
% 仍然逐候选精化；精化前后择优是本次新增的保护规则。
% 得分严格提高才采用精化直线，相等或降低时保留原直线。
x = x(:); time = time(:); ridgeTime = ridgeTime(:); ridge = ridge(:);
mList = mList(:).';
validateattributes(mList,{'numeric'},{'integer','positive','finite','nonempty'});
assert(numel(x)==numel(time) && numel(ridgeTime)==numel(ridge), ...
    '信号、时间轴或峰点的维度不一致。');
assert(all(diff(ridgeTime)>0),'峰点时间必须严格递增。');
maxM = max(mList);
if nargin<7
    pairs = wdr_sample_pairs(numel(ridge),maxM,cfg);
else
    assert(isequal(size(pairs),[2 maxM]) && all(pairs(:)==fix(pairs(:))) && ...
        all(pairs(:)>=1) && all(pairs(:)<=numel(ridge)) && ...
        all(pairs(1,:)~=pairs(2,:)),'抽样点对不合法。');
end

% 粗直线与原基线使用同样的计算，点对也可以由对照程序统一传入。
i1 = pairs(1,:); i2 = pairs(2,:);
aCoarse = ((ridge(i2)-ridge(i1))./(ridgeTime(i2)-ridgeTime(i1))).';
bCoarse = (ridge(i1)-aCoarse(:).*ridgeTime(i1)).';
aFine = zeros(1,maxM); bFine = zeros(1,maxM);
coarseScores = zeros(1,maxM); fineScores = zeros(1,maxM);
aKept = aCoarse; bKept = bCoarse;
usedRefinement = false(1,maxM);
scores = zeros(1,maxM);
for first = 1:cfg.scoreBlockSize:maxM
    block = first:min(first+cfg.scoreBlockSize-1,maxM);
    [aFine(block),bFine(block)] = wdr_oshea_refine( ...
        x,time,aCoarse(block),bCoarse(block),cfg.osheaWindowLength);
    % 两条直线都使用全部原始复样本、同一个相干得分 J 比较。
    % 此处使用非共轭转置；整体常数相位不改变相干和的幅值。
    coarsePhase = (time.^2/2)*aCoarse(block)+time*bCoarse(block);
    coarseScores(block) = abs(x.'*exp(-1j*coarsePhase));
    finePhase = (time.^2/2)*aFine(block)+time*bFine(block);
    fineScores(block) = abs(x.'*exp(-1j*finePhase));
    useFine = fineScores(block)>coarseScores(block);
    usedRefinement(block) = useFine;
    chosen = block(useFine);
    aKept(chosen) = aFine(chosen);
    bKept(chosen) = bFine(chosen);
    scores(block) = coarseScores(block);
    scores(chosen) = fineScores(chosen);
end

best = zeros(1,numel(mList));
estimate.score = zeros(1,numel(mList));
for k = 1:numel(mList)
    [estimate.score(k),best(k)] = max(scores(1:mList(k)));
end
estimate.a = aKept(best);
estimate.b = bKept(best);
estimate.omega = ridgeTime*estimate.a+ones(size(ridgeTime))*estimate.b;
estimate.mList = mList;
estimate.bestIteration = best;
estimate.bestPairs = pairs(:,best);
% 候选参数与候选得分都对应逐对择优后保留下来的直线。
estimate.candidateA = aKept;
estimate.candidateB = bKept;
estimate.candidateScores = scores;
estimate.coarseCandidateA = aCoarse;
estimate.coarseCandidateB = bCoarse;
estimate.coarseCandidateScores = coarseScores;
estimate.refinedCandidateA = aFine;
estimate.refinedCandidateB = bFine;
estimate.refinedCandidateScores = fineScores;
estimate.candidateUsedRefinement = usedRefinement;
estimate.usedRefinement = usedRefinement(best);
estimate.osheaWindowLength = cfg.osheaWindowLength;
end
