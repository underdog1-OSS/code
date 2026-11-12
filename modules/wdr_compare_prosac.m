function result = wdr_compare_prosac(cfg, growthBudget)
%WDR_COMPARE_PROSAC 在同一信号和相同 M 下比较三种抽样方式。
% 1: 原 WD-RANSAC 的高斯抽样；2: 所有时间点均匀抽样；
% 3: WD 峰点质量排序后，使用固定日程的 PROSAC 抽样。
% 不使用自适应停止或改变相干评分。所有误差均直接与真实频率相减。
% 排序质量的正确率只用于仿真诊断，不参与质量分数或估计。
validateattributes(cfg.mList, {'numeric'}, {'vector','integer','positive','finite','nonempty'});
validateattributes(cfg.trials, {'numeric'}, {'scalar','integer','positive','finite'});
validateattributes(growthBudget, {'numeric'}, {'scalar','integer','positive','finite'});
mList = cfg.mList(:).';
nM = numel(mList);
nSnr = numel(cfg.snrDb);
nTrials = cfg.trials;
result.cfg = cfg;
result.growthBudget = growthBudget;
result.snrDb = cfg.snrDb(:);
result.mList = mList;
result.gaussianError = zeros(nTrials,nSnr,nM);
result.uniformError = zeros(nTrials,nSnr,nM);
result.prosacError = zeros(nTrials,nSnr,nM);
result.wdSeconds = zeros(nTrials,nSnr);
result.qualitySeconds = zeros(nTrials,nSnr);
result.gaussianSamplingSeconds = zeros(nTrials,nSnr);
result.gaussianScoringSeconds = zeros(nTrials,nSnr);
result.uniformSamplingSeconds = zeros(nTrials,nSnr);
result.uniformScoringSeconds = zeros(nTrials,nSnr);
result.prosacSamplingSeconds = zeros(nTrials,nSnr);
result.prosacScoringSeconds = zeros(nTrials,nSnr);
result.correctAll = zeros(nTrials,nSnr);
result.correctTop10 = zeros(nTrials,nSnr);
result.correctTop20 = zeros(nTrials,nSnr);
result.trialSeed = zeros(nTrials,nSnr);
clockStart = tic;
for s = 1:nSnr
    for trial = 1:nTrials
        seed = cfg.seed+(s-1)*nTrials+trial;
        rng(seed,'twister');
        signal = wdr_generate_lfm(cfg,cfg.snrDb(s));
        timer = tic;
        [wd,omega,tWd] = wdr_distribution(signal.x,cfg);
        result.wdSeconds(trial,s) = toc(timer);
        ridge = wdr_ridge(wd,omega);
        nTime = numel(ridge);

        % 与原蒙特卡洛模块相同：信号生成后直接用高斯抽样。
        timer = tic;
        gaussianPairs = wdr_sample_pairs(nTime,max(mList),cfg);
        result.gaussianSamplingSeconds(trial,s) = toc(timer);
        timer = tic;
        gaussian = wdr_estimate(signal.x,cfg.t,tWd,ridge,mList,cfg,gaussianPairs);
        result.gaussianScoringSeconds(trial,s) = toc(timer);

        % 均匀抽样使用独立、可复现的随机序列，同一列内不能重复选点。
        rng(mod(seed+100000000,2^32),'twister');
        timer = tic;
        uniformPairs = zeros(2,max(mList));
        for k = 1:max(mList)
            uniformPairs(:,k) = sort(randperm(nTime,2));
        end
        result.uniformSamplingSeconds(trial,s) = toc(timer);
        timer = tic;
        uniform = wdr_estimate(signal.x,cfg.t,tWd,ridge,mList,cfg,uniformPairs);
        result.uniformScoringSeconds(trial,s) = toc(timer);

        rng(mod(seed+200000000,2^32),'twister');
        [prosac,quality,~,timing] = wdr_prosac_estimate( ...
            signal.x,cfg.t,tWd,ridge,wd,omega,mList,cfg,growthBudget);
        result.qualitySeconds(trial,s) = timing.qualitySeconds;
        result.prosacSamplingSeconds(trial,s) = timing.samplingSeconds;
        result.prosacScoringSeconds(trial,s) = timing.scoringSeconds;
        result.gaussianError(trial,s,:) = reshape( ...
            mean((gaussian.omega-signal.omegaTrue).^2,1),1,1,nM);
        result.uniformError(trial,s,:) = reshape( ...
            mean((uniform.omega-signal.omegaTrue).^2,1),1,1,nM);
        result.prosacError(trial,s,:) = reshape( ...
            mean((prosac.omega-signal.omegaTrue).^2,1),1,1,nM);

        % 仅在仿真中用已知真值检查排序质量；不影响抽样或模型评分。
        tolerance = 2*median(diff(omega));
        correct = abs(ridge-signal.omegaTrue)<=tolerance;
        n10 = max(1,ceil(0.10*nTime));
        n20 = max(1,ceil(0.20*nTime));
        result.correctAll(trial,s) = mean(correct);
        result.correctTop10(trial,s) = mean(correct(quality.order(1:n10)));
        result.correctTop20(trial,s) = mean(correct(quality.order(1:n20)));
        result.trialSeed(trial,s) = seed;
        if cfg.showProgress && (mod(trial,100)==0 || trial==nTrials)
            fprintf('SNR %+3g dB | %d/%d trials | elapsed %.1f s\n', ...
                cfg.snrDb(s),trial,nTrials,toc(clockStart));
        end
    end
end
% 先对线性 MSE 在试验间求平均，再转换为 dB。
result.gaussianMse = reshape(mean(result.gaussianError,1),nSnr,nM);
result.uniformMse = reshape(mean(result.uniformError,1),nSnr,nM);
result.prosacMse = reshape(mean(result.prosacError,1),nSnr,nM);
result.gaussianMseDb = 10*log10(result.gaussianMse);
result.uniformMseDb = 10*log10(result.uniformMse);
result.prosacMseDb = 10*log10(result.prosacMse);
result.correctTolerance = 2*pi/(cfg.wdLength*cfg.dt);
result.elapsedSeconds = toc(clockStart);
end
