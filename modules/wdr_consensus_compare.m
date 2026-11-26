function result = wdr_consensus_compare(cfg)
% 两方法共用信号、WD、signed 最大峰及前 B 对点；计入全部在线准备工作。
% 原 wdr_estimate 和 wdr_sample_pairs 不修改，基线仍分块评价全部 M 条直线。
if ~isfield(cfg,'consensusGlobalPeriod'), cfg.consensusGlobalPeriod=0; end
validateattributes(cfg.trials,{'numeric'},{'scalar','integer','positive','<',25000});
validateattributes(cfg.mList,{'numeric'},{'scalar','integer','positive'});
validateattributes(cfg.initialBudget,{'numeric'}, ...
    {'scalar','integer','positive','<',cfg.mList});
validateattributes(cfg.consensusTauBins,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(cfg.consensusGlobalPeriod,{'numeric'}, ...
    {'scalar','integer','nonnegative','finite'});
nSnr = numel(cfg.snrDb); nTrials = cfg.trials;
nTime = numel(cfg.centerIndices);
M = cfg.mList; B = cfg.initialBudget;
expectedExplorationCount = 0;
if cfg.consensusGlobalPeriod>0
    expectedExplorationCount = floor((M-B)/cfg.consensusGlobalPeriod);
end
result.cfg = cfg;
result.expectedGlobalExplorationCount = expectedExplorationCount;
result.snrDb = cfg.snrDb(:);
if expectedExplorationCount>0
    result.methodNames = {'原 WD-RANSAC',sprintf( ...
        '改进：后段%.0f%%全局探索＋局部反馈＋WLS',100*expectedExplorationCount/(M-B))};
else
    result.methodNames = {'原 WD-RANSAC','改进算法：反馈抽样＋局部加权＋J择优'};
end
result.error = zeros(nTrials,nSnr,2);
result.totalSeconds = zeros(nTrials,nSnr,2);
result.methodSeconds = zeros(nTrials,nSnr,2);
result.estimatedA = zeros(nTrials,nSnr,2);
result.estimatedB = zeros(nTrials,nSnr,2);
result.score = zeros(nTrials,nSnr,2); % 最终实际输出模型的 J。
result.searchBestScore = zeros(nTrials,nSnr,2);
result.searchBestA = zeros(nTrials,nSnr,2);
result.searchBestB = zeros(nTrials,nSnr,2);
result.outputScoreEvaluated = true(nTrials,nSnr,2);
result.bestIteration = zeros(nTrials,nSnr,2);
fields = {'wdSeconds','peakSeconds','gridSetupSeconds','prefixSamplingSeconds','initialBestScore', ...
    'initialA','initialB','initialError', ...
    'initialConsensusSize','finalConsensusSize','updateCount','fallbackCount', ...
    'localCandidateCount','prefixMaxScoreDifference','trueA','trueB', ...
    'globalExplorationCount','globalExplorationUpdateCount','globalCandidateCount', ...
    'localUpdateCount','fallbackUpdateCount', ...
    'trialSeed','prefixSeed','baselineSeed','improvedSeed', ...
    'weightUpdateCount','positiveWeightCount','maxWeight','weightSum','refitPositiveWeightCount', ...
    'weightedPointsOutsideFinalConsensusCount','preRefitError', ...
    'preRefitWeightedSse','postRefitWeightedSse','outputWeightedSse'};
for k = 1:numel(fields), result.(fields{k})=zeros(nTrials,nSnr); end
% 原始 WLS 候选与择优后输出分开保存；拟合退化时原始候选记为 NaN。
result.postRefitError = nan(nTrials,nSnr);
result.weightedRefitA = nan(nTrials,nSnr);
result.weightedRefitB = nan(nTrials,nSnr);
result.weightedRefitScore = nan(nTrials,nSnr);
result.candidateCount = zeros(nTrials,nSnr,2);
result.candidateScoreCount = zeros(nTrials,nSnr,2);
result.additionalScoreCount = zeros(nTrials,nSnr,2);
result.totalScoreCount = zeros(nTrials,nSnr,2);
result.initialCandidateCount = zeros(nTrials,nSnr);
result.weightedRefitApplied = false(nTrials,nSnr);
result.weightedRefitAccepted = false(nTrials,nSnr);
result.finalWeights = zeros(nTime,nTrials,nSnr,'uint32');
result.refitWeights = zeros(nTime,nTrials,nSnr,'uint32');
result.examples = cell(nSnr,1);
clockStart = tic;
for s = 1:nSnr
    for trial = 1:nTrials
        % 每次四个不重叠种子；两个方法的顺序不会改变各自抽样序列。
        dataSeed = cfg.seed+(s-1)*100000+4*trial;
        prefixSeed = dataSeed+1; baselineSeed = dataSeed+2; improvedSeed = dataSeed+3;
        rng(dataSeed,'twister');
        signal = wdr_generate_lfm(cfg,cfg.snrDb(s));
        timer = tic;
        [wd,omega,tWd] = wdr_distribution(signal.x,cfg);
        result.wdSeconds(trial,s) = toc(timer);
        timer = tic;
        ridge = wdr_ridge(wd,omega);
        result.peakSeconds(trial,s) = toc(timer);
        timer = tic;
        deltaOmega = mean(diff(omega));
        assert(all(diff(omega)>0) && max(abs(diff(omega)-deltaOmega))<1e-9, ...
            '对照需要递增均匀频率网格。');
        methodCfg = cfg;
        methodCfg.consensusTau = cfg.consensusTauBins*deltaOmega;
        result.gridSetupSeconds(trial,s) = toc(timer);
        % 前缀的随机流设置和点对生成只做一次，实际耗时完整计入两方法。
        timer = tic;
        rng(prefixSeed,'twister');
        prefixPairs = wdr_sample_pairs(numel(ridge),B,cfg);
        result.prefixSamplingSeconds(trial,s) = toc(timer);
        commonSeconds = result.wdSeconds(trial,s)+result.peakSeconds(trial,s)+result.gridSetupSeconds(trial,s)+ ...
            result.prefixSamplingSeconds(trial,s);

        order = circshift(1:2,[0 mod(s+trial,2)]);
        for method = order
            timer = tic;
            if method==1
                rng(baselineSeed,'twister');
                % 只拆开抽样准备以共用前缀；分布不变，评分仍一次处理 M 条。
                tailPairs = wdr_sample_pairs(numel(ridge),M-B,cfg);
                baselinePairs = [prefixPairs tailPairs];
                baseline = wdr_estimate(signal.x,cfg.t,tWd,ridge,M,cfg,baselinePairs);
                estimate = baseline;
            else
                rng(improvedSeed,'twister');
                improved = wdr_consensus_estimate(signal.x,cfg.t,tWd,ridge,methodCfg,prefixPairs);
                estimate = improved;
            end
            result.methodSeconds(trial,s,method) = toc(timer);
            result.totalSeconds(trial,s,method) = commonSeconds+result.methodSeconds(trial,s,method);
            % 真值误差统计在计时之外；算法模块不接收真实参数。
            result.error(trial,s,method) = mean((estimate.omega-signal.omegaTrue).^2);
            result.estimatedA(trial,s,method) = estimate.a;
            result.estimatedB(trial,s,method) = estimate.b;
            result.score(trial,s,method) = estimate.score;
            if method==1
                result.searchBestScore(trial,s,method) = estimate.score;
                result.searchBestA(trial,s,method) = estimate.a;
                result.searchBestB(trial,s,method) = estimate.b;
                result.candidateScoreCount(trial,s,method) = M;
                result.totalScoreCount(trial,s,method) = M;
            else
                result.searchBestScore(trial,s,method) = estimate.searchBestScore;
                result.searchBestA(trial,s,method) = estimate.searchBestA;
                result.searchBestB(trial,s,method) = estimate.searchBestB;
                result.outputScoreEvaluated(trial,s,method) = estimate.outputScoreEvaluated;
                result.candidateScoreCount(trial,s,method) = estimate.candidateScoreCount;
                result.additionalScoreCount(trial,s,method) = estimate.additionalScoreCount;
                result.totalScoreCount(trial,s,method) = estimate.totalScoreCount;
            end
            result.bestIteration(trial,s,method) = estimate.bestIteration;
            result.candidateCount(trial,s,method) = numel(estimate.candidateScores);
        end
        % 对比审计不计入在线估计耗时；任何断言都不按估计误差筛掉试验。
        assert(isequal(baselinePairs(:,1:B),improved.candidatePairs(:,1:B)), ...
            '两方法的初始化点对不一致。');
        prefixDifference = max(abs(baseline.candidateScores(1:B)-improved.candidateScores(1:B)));
        assert(prefixDifference<=1e-9*max(1,max(baseline.candidateScores(1:B))), ...
            '共用点对的评分不一致。');
        assert(improved.candidateCount==M && improved.initialCandidateCount==B && ...
            improved.localCandidateCount+improved.fallbackCount+ ...
            improved.globalExplorationCount==M-B && ...
            improved.globalExplorationCount==expectedExplorationCount, ...
            '实际候选数或初始化预算不正确。');
        assert(improved.weightUpdateCount==1+improved.updateCount && ...
            improved.candidateScoreCount==M && ...
            improved.additionalScoreCount==double(improved.weightedRefitApplied) && ...
            improved.totalScoreCount==M+improved.additionalScoreCount, ...
            '权重更新时机或相干评分次数不正确。');
        assert(improved.score>=improved.searchBestScore && ...
            improved.weightedRefitAccepted== ...
            (improved.weightedRefitApplied && improved.weightedRefitScore>improved.searchBestScore), ...
            '末尾 J 择优必须保留较高得分，平分时保留搜索模型。');
        result.initialCandidateCount(trial,s) = improved.initialCandidateCount;
        result.initialBestScore(trial,s) = improved.initialScore;
        result.initialA(trial,s) = improved.initialA;
        result.initialB(trial,s) = improved.initialB;
        result.initialConsensusSize(trial,s) = numel(improved.initialConsensusIndex);
        result.finalConsensusSize(trial,s) = numel(improved.finalConsensusIndex);
        result.updateCount(trial,s) = improved.updateCount;
        result.fallbackCount(trial,s) = improved.fallbackCount;
        result.localCandidateCount(trial,s) = improved.localCandidateCount;
        result.globalExplorationCount(trial,s) = improved.globalExplorationCount;
        result.globalExplorationUpdateCount(trial,s) = improved.globalExplorationUpdateCount;
        result.globalCandidateCount(trial,s) = improved.globalCandidateCount;
        result.localUpdateCount(trial,s) = improved.localUpdateCount;
        result.fallbackUpdateCount(trial,s) = improved.fallbackUpdateCount;
        result.weightUpdateCount(trial,s) = improved.weightUpdateCount;
        result.positiveWeightCount(trial,s) = nnz(improved.finalWeights);
        result.maxWeight(trial,s) = max(improved.finalWeights);
        result.weightSum(trial,s) = sum(improved.finalWeights);
        result.finalWeights(:,trial,s) = uint32(improved.finalWeights);
        result.refitWeights(:,trial,s) = uint32(improved.refitWeights);
        result.refitPositiveWeightCount(trial,s) = improved.refitPositiveWeightCount;
        result.weightedRefitApplied(trial,s) = improved.weightedRefitApplied;
        result.weightedRefitAccepted(trial,s) = improved.weightedRefitAccepted;
        result.weightedRefitA(trial,s) = improved.weightedRefitA;
        result.weightedRefitB(trial,s) = improved.weightedRefitB;
        result.weightedRefitScore(trial,s) = improved.weightedRefitScore;
        result.weightedPointsOutsideFinalConsensusCount(trial,s) = ...
            improved.weightedPointsOutsideFinalConsensusCount;
        % 初始模型是前 B 个候选中 J 最高的模型；此真值误差只用于计时外诊断。
        initialOmega = tWd*improved.initialA+improved.initialB;
        result.initialError(trial,s) = mean((initialOmega-signal.omegaTrue).^2);
        result.preRefitError(trial,s) = mean((improved.searchBestOmega-signal.omegaTrue).^2);
        if improved.weightedRefitApplied
            result.postRefitError(trial,s) = mean((improved.weightedRefitOmega-signal.omegaTrue).^2);
        end
        result.preRefitWeightedSse(trial,s) = improved.preRefitWeightedSse;
        result.postRefitWeightedSse(trial,s) = improved.postRefitWeightedSse;
        result.outputWeightedSse(trial,s) = improved.outputWeightedSse;
        result.prefixMaxScoreDifference(trial,s) = prefixDifference;
        result.trueA(trial,s) = signal.a; result.trueB(trial,s) = signal.b;
        result.trialSeed(trial,s) = dataSeed;
        result.prefixSeed(trial,s) = prefixSeed;
        result.baselineSeed(trial,s) = baselineSeed;
        result.improvedSeed(trial,s) = improvedSeed;
        if trial==1
            % 每个 SNR 保留一个完整候选历史，其余试验只留必要统计量，控制内存。
            result.examples{s} = struct('trial',trial,'prefixPairs',prefixPairs, ...
                'baselinePairs',baselinePairs,'ridge',ridge,'ridgeTime',tWd, ...
                'baseline',baseline,'improved',improved);
        end
        if cfg.showProgress && (mod(trial,100)==0 || trial==nTrials)
            fprintf('一致集反馈对照：SNR %+g dB，%d/%d 次，累计 %.1f 秒\n', ...
                cfg.snrDb(s),trial,nTrials,toc(clockStart));
        end
    end
end
result.deltaOmega = deltaOmega;
result.tau = cfg.consensusTauBins*deltaOmega;
result.omegaGrid = omega;
result.ridgeTime = tWd;
result.mse = reshape(mean(result.error,1),nSnr,2);
result.mseDb = 10*log10(result.mse);
% 与主 MSE 图相同：先在线性域平均，再转 dB；原始 WLS 退化不被悄悄删去。
result.initialMse = mean(result.initialError,1).';
result.preRefitMse = mean(result.preRefitError,1).';
result.postRefitMse = mean(result.postRefitError,1).';
result.initialMseDb = 10*log10(result.initialMse);
result.preRefitMseDb = 10*log10(result.preRefitMse);
result.postRefitMseDb = 10*log10(result.postRefitMse);
% 离线粗差诊断：轨迹 RMSE 超过一个频率格点，不反馈给任何模型选择。
% 统计次数和它们对线性域总误差的贡献，避免只看平均曲线而漏掉少数大失败。
result.grossErrorMseThreshold = deltaOmega^2;
result.grossErrorCount = reshape(sum(result.error>result.grossErrorMseThreshold,1),nSnr,2);
result.initialGrossErrorCount = sum(result.initialError>result.grossErrorMseThreshold,1).';
result.searchGrossErrorCount = sum(result.preRefitError>result.grossErrorMseThreshold,1).';
result.improvedGrossErrorMseFraction = zeros(nSnr,1);
for s = 1:nSnr
    finalErrors = result.error(:,s,2);
    totalError = sum(finalErrors);
    if totalError>0
        result.improvedGrossErrorMseFraction(s) = ...
            sum(finalErrors(finalErrors>result.grossErrorMseThreshold))/totalError;
    end
end
result.grossErrorSummary = table(result.snrDb,result.grossErrorCount(:,1), ...
    result.initialGrossErrorCount,result.searchGrossErrorCount,result.grossErrorCount(:,2), ...
    result.improvedGrossErrorMseFraction, ...
    'VariableNames',{'SNR_dB','BaselineGrossErrorCount','InitialGrossErrorCount', ...
    'SearchGrossErrorCount','OutputGrossErrorCount','GrossErrorMseFraction'});
result.meanTotalSeconds = reshape(mean(result.totalSeconds,1),nSnr,2);
result.elapsedSeconds = toc(clockStart);
result.environment.matlab = version;
result.environment.computer = computer;
result.notes = { ...
    '普通 WD-RANSAC 的高斯抽样分布与分块评分不变，只拆出共用前 B 对点。'; ...
    '共用 WD、取峰、前缀点对的实际时间均计入两组；各自剩余抽样和评分也在计时内。'; ...
    '除预定全局探索和异常回退外，后续候选按更新后的一致集逐次抽样；保持即时反馈，不提前成批生成或丢弃已评价候选。'; ...
    'consensusGlobalPeriod 为正时，每隔该次数安排一次全局探索，其他轮次正常从当前局部集抽样；为零时恢复原纯局部方案。'; ...
    '预定全局探索与不足两个一致点时的异常回退分开计数，二者都包含在剩余 M-B 个随机候选内。'; ...
    '全局探索若刷新最佳 J，同样对全部原始峰点重建一致集并更新权重；不要求新模型来自旧局部集。'; ...
    'tau 是三个实际角频率格点间隔；当前单位是 rad/归一化时间单位，不是 Hz。'; ...
    '所有 MSE 都使用中心 WD 时刻的直接角频率误差，不折回周期，不删除失败试验。'; ...
    '只记录相对当前模型的一致集；仿真真值不参与一致集、候选选择或最佳模型更新。'; ...
    '初始化最佳模型确定后记分一次，以后仅在最佳 J 刷新时记分；权重不参与随机抽样。'; ...
    '一致集始终用固定 tau；计分门限首次为 tau/3、2*tau/3、tau，后续按上次最宽加权内点集的最大残差调整。'; ...
    '未命中门限不扣分，点退出局部集也不清除其全局历史权重，但最终局部集外的点不参与末尾 WLS。'; ...
    '最终只用搜索最佳模型对应的局部集及其中各点的累计权重做一次 WLS，再用相同全信号 J 择优；严格更高才接受。'; ...
    'finalWeights 保存全部点的历史权重；refitWeights 是末尾 WLS 的有效权重，最终局部集外为零；三个阶段的加权 SSE 均使用 refitWeights。'; ...
    '随机候选仍为 M 个；有效 WLS 额外评价一次 J，拟合退化不额外评分；平分保留搜索模型。'; ...
    'searchBestScore 是搜索最佳 J，weightedRefitScore 是原始 WLS 的 J，score 是择优后输出的 J。'; ...
    '权重更新、历史记录、最终 WLS 及其 J 评分和择优均计入改进算法的总耗时。'; ...
    'initialA、initialB 和 initialError 保存前 B 次初始化结束时最佳模型的参数和逐次真值误差；initialMseDb 是先平均线性误差再转 dB 的汇总。'; ...
    'initialError、preRefitError、postRefitError 和 error(:,:,2) 分别为初始化结束、搜索结束未加权、原始 WLS、J 择优后的逐次误差，均仅在计时外用于诊断。'; ...
    '原始 WLS 拟合退化时 postRefitError 及其分数为 NaN；最终仍输出搜索模型，不删除该试验。'; ...
    'weightedRefitApplied 表示完成有效 WLS，weightedRefitAccepted 表示最终接受它；J 不下降不保证真值误差不增大。'};
result.notes{end+1} = '粗差统计仅作离线诊断：最终轨迹 RMSE 超过一个频率格点；不是系统检测失败率或捕获失败率。';
result.summary = table;
for s = 1:nSnr
    for method = 1:2
        if method==1
            meanWeightUpdates = 0; meanPositiveWeights = 0; refitFraction = 0;
            acceptedFraction = 0;
        else
            meanWeightUpdates = mean(result.weightUpdateCount(:,s));
            meanPositiveWeights = mean(result.positiveWeightCount(:,s));
            refitFraction = mean(result.weightedRefitApplied(:,s));
            acceptedFraction = mean(result.weightedRefitAccepted(:,s));
        end
        row = table(cfg.snrDb(s),string(result.methodNames{method}),nTrials,M,B, ...
            result.tau,result.mse(s,method),result.mseDb(s,method), ...
            result.meanTotalSeconds(s,method),meanWeightUpdates,meanPositiveWeights,refitFraction, ...
            acceptedFraction,mean(result.additionalScoreCount(:,s,method)), ...
            'VariableNames',{'SNR_dB','Method','Trials','M','B','Tau', ...
            'IF_MSE','IF_MSE_dB','MeanTotalSeconds', ...
            'MeanWeightUpdates','MeanPositiveWeights','WeightedRefitFraction', ...
            'WeightedRefitAcceptedFraction','MeanAdditionalScoreCount'});
        result.summary = [result.summary;row]; %#ok<AGROW>
    end
end
result.samplingSummary = table(result.snrDb, ...
    mean(result.globalExplorationCount,1).',mean(result.localCandidateCount,1).', ...
    mean(result.fallbackCount,1).',mean(result.globalExplorationUpdateCount,1).', ...
    'VariableNames',{'SNR_dB','MeanGlobalExplorationCount','MeanLocalCount', ...
    'MeanFallbackCount','MeanGlobalExplorationUpdates'});
end
