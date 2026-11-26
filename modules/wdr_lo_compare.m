function result = wdr_lo_compare(cfg)
% 师兄WD-LO-RANSAC与原WD-RANSAC的配对对照，只输出数值结构体。
% 按师兄式(4)(5)使用abs(WD)取峰，两方法共享相同峰点和全部M对样本。
% 按师兄式(1)使用复噪声总功率SNR；原信号模块通过3.0103dB换算复用。
validateattributes(cfg.trials,{'numeric'},{'scalar','integer','positive','<',25000});
validateattributes(cfg.mList,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(cfg.loTauBins,{'numeric'},{'scalar','real','positive','finite'});
validateattributes(cfg.snrDb,{'numeric'},{'vector','real','finite','nonempty'});
validateattributes(cfg.seed,{'numeric'},{'scalar','integer','nonnegative','finite'});
nTrials = cfg.trials; nSnr = numel(cfg.snrDb); M = cfg.mList;
assert(cfg.seed+(nSnr-1)*100000+4*nTrials+1<=2^32-1, ...
    '试验种子超出随机数发生器允许的范围。');
result.cfg = cfg;
result.methodNames = {'WD-RANSAC基线','师兄 WD-LO-RANSAC'};
result.snrDb = cfg.snrDb(:);
result.snrConvention = '信号功率/复噪声总功率';
result.peakConvention = '每个时刻取abs(WD)的最大峰，两方法相同';
result.componentSnrDb = result.snrDb+10*log10(2);
result.error = zeros(nTrials,nSnr,2);
result.totalSeconds = zeros(nTrials,nSnr,2);
result.methodSeconds = zeros(nTrials,nSnr,2);
result.estimatedA = zeros(nTrials,nSnr,2);
result.estimatedB = zeros(nTrials,nSnr,2);
result.score = zeros(nTrials,nSnr,2);
result.candidateCount = M*ones(nTrials,nSnr,2);
result.additionalScoreCount = zeros(nTrials,nSnr,2);
result.loTriggerCount = zeros(nTrials,nSnr);
result.loFitCount = zeros(nTrials,nSnr);
result.loAcceptedCount = zeros(nTrials,nSnr);
result.outputFromLo = false(nTrials,nSnr);
result.wdSeconds = zeros(nTrials,nSnr);
result.peakSeconds = zeros(nTrials,nSnr);
result.samplingSeconds = zeros(nTrials,nSnr);
result.gridSetupSeconds = zeros(nTrials,nSnr);
result.trueA = zeros(nTrials,nSnr); result.trueB = zeros(nTrials,nSnr);
result.trialSeed = zeros(nTrials,nSnr); result.samplingSeed = zeros(nTrials,nSnr);
result.orderFirstMethod = zeros(nTrials,nSnr);
result.candidateMaxScoreDifference = zeros(nTrials,nSnr);
result.examples = cell(nSnr,1);
clockStart = tic;

for s = 1:nSnr
    for trial = 1:nTrials
        dataSeed = cfg.seed+(s-1)*100000+4*trial;
        samplingSeed = dataSeed+1;
        rng(dataSeed,'twister');
        % 原模块每个实分量方差为A^2*10^(-SNR/10)。加3.0103dB后，
        % 其两个实分量方差之和恰为本次论文SNR要求的复噪声总功率。
        signal = wdr_generate_lfm(cfg,result.componentSnrDb(s));

        timer = tic;
        [wd,omega,tWd] = wdr_distribution(signal.x,cfg);
        result.wdSeconds(trial,s) = toc(timer);
        timer = tic;
        [ridge,amplitude] = wdr_ridge(abs(wd),omega);
        result.peakSeconds(trial,s) = toc(timer);
        timer = tic;
        deltaOmega = mean(diff(omega));
        assert(all(diff(omega)>0) && max(abs(diff(omega)-deltaOmega))<1e-9, ...
            '需要递增且等间隔的角频率网格。');
        methodCfg = cfg; methodCfg.loTau = cfg.loTauBins*deltaOmega;
        result.gridSetupSeconds(trial,s) = toc(timer);
        timer = tic;
        rng(samplingSeed,'twister');
        pairs = wdr_sample_pairs(numel(ridge),M,cfg);
        result.samplingSeconds(trial,s) = toc(timer);
        % 共用操作只实际计算一次，但其耗时完整计入两种算法各自的总时间。
        commonSeconds = result.wdSeconds(trial,s)+result.peakSeconds(trial,s)+ ...
            result.gridSetupSeconds(trial,s)+result.samplingSeconds(trial,s);
        order = circshift(1:2,[0 mod(s+trial,2)]);
        result.orderFirstMethod(trial,s) = order(1);
        for method = order
            timer = tic;
            if method==1
                baseline = wdr_estimate(signal.x,cfg.t,tWd,ridge,M,cfg,pairs);
                estimate = baseline;
            else
                improved = wdr_lo_estimate( ...
                    signal.x,cfg.t,tWd,ridge,amplitude,methodCfg,pairs);
                estimate = improved;
            end
            result.methodSeconds(trial,s,method) = toc(timer);
            result.totalSeconds(trial,s,method) = commonSeconds+result.methodSeconds(trial,s,method);
            % 真值仅用于计时外评价；保留全部失败试验，使用直接误差。
            result.error(trial,s,method) = mean((estimate.omega-signal.omegaTrue).^2);
            result.estimatedA(trial,s,method) = estimate.a;
            result.estimatedB(trial,s,method) = estimate.b;
            result.score(trial,s,method) = estimate.score;
        end
        % 一致性审计位于计时之外，不用于挑选输出或删去失败试验。
        scoreDifference = max(abs(baseline.candidateScores-improved.candidateScores));
        assert(scoreDifference<=1e-9*max(1,max(baseline.candidateScores)), ...
            '同一批原始候选的J计算不一致。');
        assert(improved.triggerFlags(1) && improved.candidateCount==M && ...
            improved.loFitCount<=improved.loTriggerCount && ...
            improved.totalScoreCount==M+improved.loFitCount, ...
            '局部优化触发或实际评分次数不正确。');
        assert(all(diff(improved.bestScoreHistory)>=0) && improved.score>=baseline.score, ...
            '必须保留所有已评价原始候选和重估候选中的最高J。');
        result.candidateMaxScoreDifference(trial,s) = scoreDifference;
        result.loTriggerCount(trial,s) = improved.loTriggerCount;
        result.loFitCount(trial,s) = improved.loFitCount;
        result.loAcceptedCount(trial,s) = improved.loAcceptedCount;
        result.outputFromLo(trial,s) = improved.outputFromLo;
        result.additionalScoreCount(trial,s,2) = improved.additionalScoreCount;
        result.trueA(trial,s) = signal.a; result.trueB(trial,s) = signal.b;
        result.trialSeed(trial,s) = dataSeed;
        result.samplingSeed(trial,s) = samplingSeed;
        if trial==1
            result.examples{s} = struct('trial',trial,'ridgeTime',tWd,'ridge',ridge, ...
                'amplitude',amplitude,'pairs',pairs,'baseline',baseline,'improved',improved);
        end
        if cfg.showProgress && (mod(trial,100)==0 || trial==nTrials)
            fprintf('师兄算法对照：SNR %+g dB，%d/%d次，累计%.1f秒\n', ...
                cfg.snrDb(s),trial,nTrials,toc(clockStart));
        end
    end
end
result.omegaGrid = omega; result.ridgeTime = tWd;
result.deltaOmega = deltaOmega; result.tau = cfg.loTauBins*deltaOmega;
result.mse = reshape(mean(result.error,1),nSnr,2);
result.mseDb = 10*log10(result.mse);
result.meanTotalSeconds = reshape(mean(result.totalSeconds,1),nSnr,2);
result.meanLoTriggers = mean(result.loTriggerCount,1).';
result.meanLoFits = mean(result.loFitCount,1).';
result.elapsedSeconds = toc(clockStart);
result.environment = struct('matlab',version,'computer',computer);
result.summary = table(result.snrDb,result.mseDb(:,1),result.mseDb(:,2), ...
    result.meanTotalSeconds(:,1),result.meanTotalSeconds(:,2), ...
    result.meanLoTriggers,result.meanLoFits, ...
    'VariableNames',{'SNR_dB','BaselineMSE_dB','LoMSE_dB', ...
    'BaselineSeconds','LoSeconds','MeanLoTriggers','MeanExtraJScores'});
result.notes = { ...
    '依据师兄稿件第2.3节式(8)-(10)：历史最优触发，一致点幅度与残差双权重，一次WLS。'; ...
    '每次重新计算权重，不累计权重；全部随机候选都从全局峰点产生，不做局部反馈抽样。'; ...
    '两方法共享全部M对点及相同WD峰点，原始候选均调用原wdr_estimate的分块评分。'; ...
    '师兄式(4)(5)使用abs(WD)取峰；旧run_demo和run_monte_carlo仍保留原signed取峰接口。'; ...
    '本入口SNR使用复噪声总功率定义，比旧模块实部或虚部分量方差定义低3.0103dB。'; ...
    '稿件判决基准存在歧义：本实现先保留原候选，再仅接受更高J的WLS结果，平分保留原候选。'; ...
    'WLS采用时间中心化和缩放后的sqrt(w)线性求解，与稿件加权平方误差目标等价。'; ...
    '退化情况包括有效正权重点不足、时间秩不足和零幅度；均保留原候选，不额外重抽。'; ...
    '每个有效局部重估额外评价一次全信号J，不把额外J计算藏在M个随机候选之内。'; ...
    '抽样沿用原高斯抽样；真参数仍按允许区域均匀生成，与稿件矩形均匀抽样后拒绝越界的分布等价。'; ...
    '现有WD的偶数长度延迟端点采用共轭平均；稿件未明确偶数端点数值细节，因此不能保证逐点重现原图。'; ...
    '共享WD、取峰、网格设置和点对生成耗时分别计入两方法，LO方法另含扫描、权重、WLS及额外J。'; ...
    '直接角频率误差在时间与试验上线性平均后转dB，不作周期折回，不以真值筛选候选。'; ...
    '采样间隔按稿件表1解释为秒，角频率为rad/s，a为rad/s^2，b为t=0时的角频偏。'; ...
    '稿件使用10000次试验，本入口默认1000次；论文的精度和耗时收益需实际运行后验证。'};
end
