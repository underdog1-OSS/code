function estimate = wdr_lo_estimate(x,time,ridgeTime,ridge,amplitude,cfg,pairs)
% 师兄稿件第2.3节：历史最优触发、一致点双权重、一次WLS重估。
% 所有随机候选仍来自全局峰点，局部优化不改变后续抽样来源。
% 原始候选继续调用wdr_estimate分块评分；仅扫描得分较高的候选做额外处理。
% 每次权重由当前候选的幅度和残差重新计算，不使用历史累计权重。
% 原文图1的判决基准有歧义：本实现先保留刷新纪录的原候选，
% 再仅在重估候选J严格更高时替换，确保保留已评价模型中的最高J。
x = x(:); time = time(:); ridgeTime = ridgeTime(:);
ridge = ridge(:); amplitude = amplitude(:);
M = cfg.mList; tau = cfg.loTau;
validateattributes(M,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(tau,{'numeric'},{'scalar','real','positive','finite'});
validateattributes(amplitude,{'numeric'},{'real','finite','nonnegative','nonempty'});
assert(numel(amplitude)==numel(ridge) && numel(ridge)==numel(ridgeTime), ...
    '峰点频率、时间和幅度的长度必须相同。');
assert(numel(ridge)>=2 && all(diff(ridgeTime)>0), ...
    '需要至少两个不同时刻的峰点。');
if nargin<7
    pairs = wdr_sample_pairs(numel(ridge),M,cfg);
end
raw = wdr_estimate(x,time,ridgeTime,ridge,M,cfg,pairs);
bestScore = -Inf; bestA = NaN; bestB = NaN;
bestIteration = 0; outputFromLo = false;
triggerFlags = false(1,M); fitFlags = false(1,M); acceptedFlags = false(1,M);
loScores = nan(1,M); loA = nan(1,M); loB = nan(1,M);
inlierCount = zeros(1,M); positiveWeightCount = zeros(1,M);
bestScoreHistory = zeros(1,M);
halfTimeSquared = time.^2/2;

for m = 1:M
    if raw.candidateScores(m)>bestScore
        triggerFlags(m) = true;
        % 先保存本轮原始候选，防止WLS退步时丢失已找到的更高J模型。
        bestA = raw.candidateA(m); bestB = raw.candidateB(m);
        bestScore = raw.candidateScores(m);
        bestIteration = m; outputFromLo = false;
        residual = abs(ridge-(bestA*ridgeTime+bestB));
        members = find(residual<=tau);
        inlierCount(m) = numel(members);
        if numel(members)>=2
            maxAmplitude = max(amplitude(members));
            if maxAmplitude>0
                % 式(9)：当前一致集内幅度归一化，乘上残差的双平方权重。
                u = residual(members)/tau;
                weights = (amplitude(members)/maxAmplitude).*(1-u.^2).^2;
                [aLo,bLo,fit] = wdr_weighted_line_fit( ...
                    ridgeTime(members),ridge(members),weights);
                positiveWeightCount(m) = fit.positivePointCount;
                if fit.applied
                    fitFlags(m) = true;
                    loA(m) = aLo; loB(m) = bLo;
                    % 与原候选一致，使用全部原始信号样本计算式(7)的J。
                    loScores(m) = abs(x.'*exp(-1j*(halfTimeSquared*aLo+time*bLo)));
                    if loScores(m)>bestScore
                        bestA = aLo; bestB = bLo; bestScore = loScores(m);
                        acceptedFlags(m) = true; outputFromLo = true;
                    end
                end
            end
        end
    end
    bestScoreHistory(m) = bestScore;
end

% 保留原始候选历史，便于检查两方法是否使用了同一批随机候选。
estimate = raw;
estimate.a = bestA; estimate.b = bestB; estimate.score = bestScore;
estimate.omega = ridgeTime*bestA+bestB;
estimate.bestIteration = bestIteration;
estimate.bestPairs = pairs(:,bestIteration);
estimate.outputFromLo = outputFromLo;
estimate.rawBestScore = raw.score;
estimate.rawBestA = raw.a; estimate.rawBestB = raw.b;
estimate.rawBestIteration = raw.bestIteration;
estimate.triggerFlags = triggerFlags;
estimate.fitFlags = fitFlags;
estimate.acceptedFlags = acceptedFlags;
estimate.loCandidateScores = loScores;
estimate.loCandidateA = loA; estimate.loCandidateB = loB;
estimate.inlierCount = inlierCount;
estimate.positiveWeightCount = positiveWeightCount;
estimate.bestScoreHistory = bestScoreHistory;
estimate.loTriggerCount = nnz(triggerFlags);
estimate.loFitCount = nnz(fitFlags);
estimate.loAcceptedCount = nnz(acceptedFlags);
estimate.candidateCount = M;
estimate.additionalScoreCount = estimate.loFitCount;
estimate.totalScoreCount = M+estimate.additionalScoreCount;
estimate.tau = tau;
end
