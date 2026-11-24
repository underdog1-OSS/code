function result = wdr_compare_oshea(cfg)
% 在相同数据和抽样预算下比较原 WD-RANSAC 与逐候选 O'Shea 精化版本。
% 两组均取全频段最大峰；精化组逐对比较原直线和精化直线的相干得分。
validateattributes(cfg.trials,{'numeric'},{'scalar','integer','positive','<',100000});
validateattributes(cfg.mList,{'numeric'},{'scalar','integer','positive'});
validateattributes(cfg.osheaWindowLength,{'numeric'}, ...
    {'scalar','integer','positive','<=',cfg.nSamples-2});
assert(mod(cfg.osheaWindowLength,2)==1,'移动平均长度必须是奇数。');
nTrials = cfg.trials; nSnr = numel(cfg.snrDb);
result.cfg = cfg;
result.snrDb = cfg.snrDb(:);
result.methodNames = {'原 WD-RANSAC','WD-RANSAC + O''Shea 评分择优'};
result.error = zeros(nTrials,nSnr,2);
result.estimatedA = zeros(nTrials,nSnr,2);
result.estimatedB = zeros(nTrials,nSnr,2);
result.score = zeros(nTrials,nSnr,2);
result.bestIteration = zeros(nTrials,nSnr,2);
result.estimationSeconds = zeros(nTrials,nSnr,2);
result.totalSeconds = zeros(nTrials,nSnr,2);
result.wdSeconds = zeros(nTrials,nSnr);
result.ridgeSeconds = zeros(nTrials,nSnr);
result.pairSeconds = zeros(nTrials,nSnr);
result.trueA = zeros(nTrials,nSnr);
result.trueB = zeros(nTrials,nSnr);
result.trialSeed = zeros(nTrials,nSnr);
clockStart = tic;
for s = 1:nSnr
    for trial = 1:nTrials
        seed = cfg.seed+(s-1)*100000+trial;
        rng(seed,'twister');
        % 沿用基线在允许参数区域内随机生成 a、b 的原始信号模型。
        signal = wdr_generate_lfm(cfg,cfg.snrDb(s));
        timer = tic;
        [wd,omega,tWd] = wdr_distribution(signal.x,cfg);
        result.wdSeconds(trial,s) = toc(timer);
        timer = tic;
        ridge = wdr_ridge(wd,omega);
        result.ridgeSeconds(trial,s) = toc(timer);
        timer = tic;
        pairs = wdr_sample_pairs(numel(ridge),cfg.mList,cfg);
        result.pairSeconds(trial,s) = toc(timer);
        result.trueA(trial,s) = signal.a;
        result.trueB(trial,s) = signal.b;
        result.trialSeed(trial,s) = seed;
        commonSeconds = result.wdSeconds(trial,s)+ ...
            result.ridgeSeconds(trial,s)+result.pairSeconds(trial,s);

        % 交替两组执行顺序，避免始终让某一组先运行。
        order = circshift(1:2,[0 mod(trial+s,2)]);
        for m = order
            timer = tic;
            if m==1
                estimate = wdr_estimate(signal.x,cfg.t,tWd,ridge,cfg.mList,cfg,pairs);
            else
                estimate = wdr_oshea_estimate(signal.x,cfg.t,tWd,ridge,cfg.mList,cfg,pairs);
            end
            result.estimationSeconds(trial,s,m) = toc(timer);
            % 两组总时间都包含 WD、全频取峰、抽样、拟合及评分；
            % 精化组额外包含全部候选的去啁啾、滤波、展开、重估和前后评分比较。
            result.totalSeconds(trial,s,m) = commonSeconds+result.estimationSeconds(trial,s,m);
            result.error(trial,s,m) = mean((estimate.omega-signal.omegaTrue).^2);
            result.estimatedA(trial,s,m) = estimate.a;
            result.estimatedB(trial,s,m) = estimate.b;
            result.score(trial,s,m) = estimate.score;
            result.bestIteration(trial,s,m) = estimate.bestIteration;
        end
        % 失效试验也参与平均，不根据真值、评分或误差删除任何试验。
        if cfg.showProgress && (mod(trial,20)==0 || trial==nTrials)
            fprintf('O''Shea 对照：SNR %+g dB，%d/%d 次，累计 %.1f 秒\n', ...
                cfg.snrDb(s),trial,nTrials,toc(clockStart));
        end
    end
end
% MSE 先在直接角频率平方误差上平均，再转 dB，与原基线口径一致。
result.mse = reshape(mean(result.error,1),nSnr,2);
result.mseDb = 10*log10(result.mse);
result.meanTotalMilliseconds = 1000*reshape(mean(result.totalSeconds,1),nSnr,2);
result.elapsedSeconds = toc(clockStart);
result.summary = table;
for s = 1:nSnr
    for m = 1:2
        row = table(cfg.snrDb(s),string(result.methodNames{m}),cfg.mList,cfg.trials, ...
            result.mse(s,m),result.mseDb(s,m),result.meanTotalMilliseconds(s,m), ...
            'VariableNames',{'SNR_dB','Method','M','Trials','IF_MSE','IF_MSE_dB','MeanTotal_ms'});
        result.summary = [result.summary;row]; %#ok<AGROW>
    end
end
end
