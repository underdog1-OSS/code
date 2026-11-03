function result = wdr_monte_carlo(cfg)
% 实验模块：在各 SNR 下重复生成随机信号、估计频率，并汇总 MSE。
%WDR_MONTE_CARLO Example 1, random admissible LFM parameters in every trial.
% Means are formed in linear error units before conversion to dB.
% Same noisy signal/ridge and nested candidate sequence for all M per trial.
% Per-trial RNG seeds make any particular failure reproducible independently.
nSnr = numel(cfg.snrDb);
nM = numel(cfg.mList);
result.cfg = cfg;
result.snrDb = cfg.snrDb(:);
result.wdError = zeros(cfg.trials,nSnr);                                    % 某次试验、某个 SNR 下，WD 直接提取频率的误差
result.ransacError = zeros(cfg.trials,nSnr,nM);                             % 同一次试验和SNR下，第m种M设置的WDRANSAC误差
result.trueA = zeros(cfg.trials,nSnr);                                      % 生成信号时使用的真实斜率 a 和截距 b
result.trueB = zeros(cfg.trials,nSnr);                                     
result.estimatedA = zeros(cfg.trials,nSnr,nM);                              % WDRANSAC 在各个 M 设置下估计出的 a 和 b
result.estimatedB = zeros(cfg.trials,nSnr,nM);
result.trialSeed = zeros(cfg.trials,nSnr);                                  % 随机种子、WD 计算时间、WDRANSAC 计算时间
result.wdSeconds = zeros(cfg.trials,nSnr);
result.ransacSeconds = zeros(cfg.trials,nSnr); % All M budgets together
clockStart = tic;                                                           % 启动整个实验的计时器
for s = 1:nSnr                                                              % 先选一个 SNR，再在这个 SNR 下重复 cfg.trials 次试验            
    for trial = 1:cfg.trials
        seed = cfg.seed + (s-1)*cfg.trials + trial;                         % 第一行是给当前这次试验编一个不同的“随机种子”
        rng(seed, 'twister');                                               % 决定后续 rand、randn 会产生哪一串随机数，相同种子生成相同随机数
        signal = wdr_generate_lfm(cfg,cfg.snrDb(s));
        timer = tic;
        [wd,omega,tWd] = wdr_distribution(signal.x,cfg);
        result.wdSeconds(trial,s) = toc(timer);
        ridge = wdr_ridge(wd,omega);
        timer = tic;
        estimate = wdr_estimate(signal.x,cfg.t,tWd,ridge,cfg.mList,cfg);
        result.ransacSeconds(trial,s) = toc(timer);
        % 原来的 WD 误差计算（直接相减）；恢复时取消本行注释并注释下面四行。
        % result.wdError(trial,s) = mean((ridge-signal.omegaTrue).^2);
        % 对照实验：把 WD 频率误差折回一个周期，观察频率边界对 MSE 的影响。
        wdPeriod = pi/cfg.dt;
        wdDifference = ridge-signal.omegaTrue;
        wdDifference = mod(wdDifference+wdPeriod/2,wdPeriod)-wdPeriod/2;
        result.wdError(trial,s) = mean(wdDifference.^2);
        % WDRANSAC 使用直接频率误差，避免把严重偏离的直线折回。
        result.ransacError(trial,s,:) = reshape( ...
            mean((estimate.omega-signal.omegaTrue).^2,1),1,1,nM);
        % 周期误差的对照版本，保留供查看；启用时需注释上面两行。
        % ransacPeriod = pi/cfg.dt;
        % ransacDifference = estimate.omega-signal.omegaTrue;
        % ransacDifference = mod(ransacDifference+ransacPeriod/2,ransacPeriod)-ransacPeriod/2;
        % result.ransacError(trial,s,:) = reshape( ...
        %     mean(ransacDifference.^2,1),1,1,nM);
        result.trueA(trial,s) = signal.a;                                   % 把这次试验的参数和种子存起来：
        result.trueB(trial,s) = signal.b;
        result.estimatedA(trial,s,:) = reshape(estimate.a,1,1,nM);
        result.estimatedB(trial,s,:) = reshape(estimate.b,1,1,nM);
        result.trialSeed(trial,s) = seed;
        if cfg.showProgress && (mod(trial,100)==0 || trial==cfg.trials)     % 每100次实验打印一次进度
            fprintf('SNR %+3g dB | %d/%d trials | elapsed %.1f s\n', ...
                cfg.snrDb(s),trial,cfg.trials,toc(clockStart));
        end
    end
end
result.wdMse = mean(result.wdError,1).';                                    % 对第一列即固定snr的所有次实验误差求平均
result.ransacMse = reshape(mean(result.ransacError,1),nSnr,nM);             % 每个SNR与M的组合得到一个平均误差，最终排成nSnr×nM矩阵；行对应SNR，列对应M
result.wdMseDb = 10*log10(result.wdMse);                                    % 把平均后的 MSE换成 dB
result.ransacMseDb = 10*log10(result.ransacMse);
result.elapsedSeconds = toc(clockStart);                                    % 记录整个蒙特卡洛实验的总用时
end

