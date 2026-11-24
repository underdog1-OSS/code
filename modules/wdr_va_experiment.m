function result = wdr_va_experiment(cfg)
% 两层独立验证和训练测试；所有路径都由 WD＋F^RR 基本 VA 得到。
result.cfg = cfg;
result.notes = cfg.notes;
result.environment.matlab = version;
result.environment.computer = computer;
result.environment.processor = getenv('PROCESSOR_IDENTIFIER');
clockStart = tic;

% 第一层：论文明确给出的线性相位，时频表示按本次要求用现有 WD。
fprintf('\n第一层：论文线性相位＋现有 WD。训练 lambda，目标为原始脊线 MSE。\n');
result.paperTraining = train_lambda(cfg,cfg.paperParameters, ...
    cfg.paperTrainSnrDb,cfg.paperTrainTrials,cfg.seed+3000000,'raw');
fprintf('第一层固定 lambda = %.7g。\n',result.paperTraining.lambda);
result.paper = run_case(cfg,cfg.paperParameters,cfg.paperSnrDb,cfg.paperTrials, ...
    cfg.seed+4000000,result.paperTraining.lambda,'论文线性信号＋现有 WD');

% 第二层：原基线的随机允许参数分布；仅使用独立训练样本挑参数。
fprintf('\n第二层：原 WD-RANSAC 随机信号模型。训练目标为直线拟合后 MSE。\n');
result.training = train_lambda(cfg,[],cfg.vaTrainSnrDb,cfg.vaTrainTrials, ...
    cfg.seed+1000000,'fit');
fprintf('统一对照固定 lambda = %.7g；测试阶段不再选择参数。\n',result.training.lambda);
% 短试跑提供当前机器的耗时证据，不参与正式统计，不改变固定参数。
result.pilot = run_case(cfg,[],[-5 20],cfg.pilotTrials,cfg.seed+8000000, ...
    result.training.lambda,'正式统计前短试跑');
fprintf('试跑单方法平均总时间：WD-RANSAC %.4f 秒，WD＋VA＋OLS %.4f 秒。\n', ...
    mean(result.pilot.totalSeconds(:,:,1),'all'),mean(result.pilot.totalSeconds(:,:,2),'all'));
result.test = run_case(cfg,[],cfg.snrDb,cfg.trials,cfg.seed+2000000, ...
    result.training.lambda,'独立统一 WD 测试');
result.elapsedSeconds = toc(clockStart);
end

function training = train_lambda(cfg,parameters,snrDb,trials,seedBase,target)
% 同一个训练样本在各 lambda 间复用，不与测试样本混用。
grid = cfg.vaLambdaGrid(:).';
rawMse = zeros(trials,numel(snrDb),numel(grid));
fitMse = rawMse; seeds = zeros(trials,numel(snrDb));
for s = 1:numel(snrDb)
    for trial = 1:trials
        seed = seedBase+s*10000+trial; rng(seed,'twister');
        signal = wdr_generate_lfm(cfg,snrDb(s),parameters);
        [wd,omega,time] = wdr_distribution(signal.x,cfg);
        truth = signal.a*time+signal.b;
        design = [time ones(size(time))];
        for k = 1:numel(grid)
            va = wdr_va_ridge(wd,omega,grid(k),cfg.vaFloorRelative);
            coefficients = design\va.omega;
            rawMse(trial,s,k) = mean((va.omega-truth).^2);
            fitMse(trial,s,k) = mean((design*coefficients-truth).^2);
        end
        seeds(trial,s) = seed;
    end
    fprintf('lambda 训练：SNR %+g dB，%d 次，%d 个预定候选。\n',snrDb(s),trials,numel(grid));
end
if strcmp(target,'raw'), loss = rawMse; else, loss = fitMse; end
objective = reshape(mean(mean(loss,1),2),1,numel(grid));
[~,best] = min(objective);
training.lambda = grid(best);
training.grid = grid;
training.objective = objective;
training.target = target;
training.rawMse = rawMse;
training.fitMse = fitMse;
training.trialSeed = seeds;
training.snrDb = snrDb;
training.trials = trials;
training.atGridBoundary = best==1 || best==numel(grid);
training.selection = '仅训练数据的线性域 MSE 在各 SNR 等权平均，选一个固定 lambda；测试不调参。';
end

function r = run_case(cfg,parameters,snrDb,trials,seedBase,lambda,label)
% 原 WD-RANSAC 和 WD＋VA＋OLS 配对比较，每次只生成一份含噪信号和 WD。
nSnr = numel(snrDb); nTime = numel(cfg.centerIndices);
r.label = label; r.lambda = lambda; r.snrDb = snrDb(:); r.trials = trials;
r.methodNames = {'原 WD-RANSAC','WD＋VA 全局脊线＋OLS'};
r.msePerTrial = zeros(trials,nSnr,2);
r.aHat = zeros(trials,nSnr,2); r.bHat = r.aHat;
r.totalSeconds = r.aHat;
r.rawVaMse = zeros(trials,nSnr);
r.peakMse = r.rawVaMse;
r.gridOracleRawMse = r.rawVaMse; r.gridOracleFitMse = r.rawVaMse;
r.vaCorrectFraction = r.rawVaMse; r.vaMeanAbsJumpBins = r.rawVaMse;
r.truthNearEdgeFraction = r.rawVaMse;
r.trueA = r.rawVaMse; r.trueB = r.rawVaMse;
r.trialSeed = r.rawVaMse; r.wdSeconds = r.rawVaMse;
r.pathIndex = zeros(nTime,trials,nSnr,'uint16');
r.examples = cell(nSnr,1);
for s = 1:nSnr
    for trial = 1:trials
        seed = seedBase+s*10000+trial; rng(seed,'twister');
        signal = wdr_generate_lfm(cfg,snrDb(s),parameters);
        timer = tic;
        [wd,omega,time] = wdr_distribution(signal.x,cfg);
        wdSeconds = toc(timer);
        truth = signal.a*time+signal.b;
        % 两组执行顺序交替；各自的全部估计工作都放在对应计时器中。
        order = circshift(1:2,[0 mod(s+trial,2)]);
        for method = order
            timer = tic;
            if method==1
                peak = wdr_ridge(wd,omega);
                estimate = wdr_estimate(signal.x,cfg.t,time,peak,cfg.mList,cfg);
                aHat = estimate.a; bHat = estimate.b; fit = estimate.omega;
            else
                va = wdr_va_ridge(wd,omega,lambda,cfg.vaFloorRelative);
                design = [time ones(size(time))];
                coefficients = design\va.omega;
                aHat = coefficients(1); bHat = coefficients(2);
                fit = design*coefficients;
            end
            r.totalSeconds(trial,s,method) = wdSeconds+toc(timer);
            r.aHat(trial,s,method) = aHat;
            r.bHat(trial,s,method) = bHat;
            r.msePerTrial(trial,s,method) = mean((fit-truth).^2);
        end
        % 下面的仿真真值仅用于离线诊断，绝不反馈给 VA 或测试参数选择。
        delta = mean(diff(omega));
        oracleIndex = round((truth-omega(1))/delta)+1;
        oracleIndex = min(numel(omega),max(1,oracleIndex));
        oracle = omega(oracleIndex);
        design = [time ones(size(time))];
        oracleFit = design*(design\oracle);
        r.rawVaMse(trial,s) = mean((va.omega-truth).^2);
        r.peakMse(trial,s) = mean((peak-truth).^2);
        r.gridOracleRawMse(trial,s) = mean((oracle-truth).^2);
        r.gridOracleFitMse(trial,s) = mean((oracleFit-truth).^2);
        r.vaCorrectFraction(trial,s) = mean(abs(va.omega-truth)<=delta);
        r.vaMeanAbsJumpBins(trial,s) = mean(abs(diff(va.index)));
        r.truthNearEdgeFraction(trial,s) = mean(truth<omega(1)+delta | truth>omega(end)-delta);
        r.trueA(trial,s) = signal.a; r.trueB(trial,s) = signal.b;
        r.trialSeed(trial,s) = seed; r.wdSeconds(trial,s) = wdSeconds;
        r.pathIndex(:,trial,s) = uint16(va.index);
        % 真 IF 为直线时，同时间区间内的正交投影不应增大轨迹 MSE。
        assert(r.msePerTrial(trial,s,2)<=r.rawVaMse(trial,s)+1e-8*max(1,r.rawVaMse(trial,s)), ...
            '直线拟合、真值或评价时间轴不一致。');
        if trial==1
            r.examples{s} = struct('seed',seed,'x',signal.x,'wd',wd,'omega',omega, ...
                'time',time,'truth',truth,'peak',peak,'va',va);
        end
        if cfg.showProgress && (mod(trial,20)==0 || trial==trials)
            fprintf('%s：SNR %+g dB，%d/%d 次。\n',label,snrDb(s),trial,trials);
        end
    end
end
r.omega = omega; r.time = time;
r.deltaOmega = mean(diff(omega));
r.deltaHz = r.deltaOmega/(2*pi*cfg.secondsPerTimeUnit);
r.lambdaOmega = lambda/r.deltaOmega^2;
r.lambdaHz = lambda/r.deltaHz^2;
r.mse = reshape(mean(r.msePerTrial,1),nSnr,2);
r.mseDb = 10*log10(r.mse);
r.mseHz2 = r.mse/(2*pi*cfg.secondsPerTimeUnit)^2;
r.aHatHzPerSec = r.aHat/(2*pi*cfg.secondsPerTimeUnit^2);
r.bHatHz = r.bHat/(2*pi*cfg.secondsPerTimeUnit);
r.frequencyAtFirstSampleHz = (r.bHat+cfg.t(1)*r.aHat)/(2*pi*cfg.secondsPerTimeUnit);
r.meanTotalSeconds = reshape(mean(r.totalSeconds,1),nSnr,2);
r.successFraction = reshape(mean(r.msePerTrial<=(cfg.successRmseBins*r.deltaOmega)^2,1),nSnr,2);
r.thresholdBracketDb = nan(2,2);
for m = 1:2
    for s = 1:nSnr
        if all(r.successFraction(s:end,m)>=cfg.successRateTarget)
            if s>1, r.thresholdBracketDb(m,1)=snrDb(s-1); else, r.thresholdBracketDb(m,1)=-Inf; end
            r.thresholdBracketDb(m,2)=snrDb(s);
            break;
        end
    end
end
r.summary = table;
for s = 1:nSnr
    for m = 1:2
        aError = r.aHat(:,s,m)-r.trueA(:,s);
        bError = r.bHat(:,s,m)-r.trueB(:,s);
        row = table(snrDb(s),string(r.methodNames{m}),trials,lambda,r.mse(s,m),r.mseDb(s,m), ...
            r.mseHz2(s,m),r.meanTotalSeconds(s,m),r.successFraction(s,m), ...
            mean(aError),sqrt(mean(aError.^2)),mean(bError),sqrt(mean(bError.^2)), ...
            'VariableNames',{'SNR_dB','Method','Trials','LambdaBin','MSE_Angular','MSE_dB', ...
            'MSE_Hz2','MeanTotalSeconds','SuccessFraction','A_Bias','A_RMSE','B_Bias','B_RMSE'});
        r.summary = [r.summary;row]; %#ok<AGROW>
    end
end
r.diagnostics = table(snrDb(:),mean(r.rawVaMse,1).',r.mse(:,2), ...
    mean(r.peakMse,1).',mean(r.gridOracleRawMse,1).',mean(r.gridOracleFitMse,1).', ...
    mean(r.vaCorrectFraction,1).',mean(r.vaMeanAbsJumpBins,1).',mean(r.truthNearEdgeFraction,1).', ...
    'VariableNames',{'SNR_dB','RawVaMSE','FitVaMSE','MaxPeakMSE','OracleGridRawMSE', ...
    'OracleGridFitMSE','VaCorrectFraction','MeanAbsJumpBins','TruthNearEdgeFraction'});
% 以整次试验为重抽样单位，报告配对 MSE 差的不确定性。
state = rng; cleaner = onCleanup(@() rng(state)); %#ok<NASGU>
rng(seedBase+9000,'twister');
r.pairedMseDifferenceCI = zeros(nSnr,3);
for s = 1:nSnr
    difference = r.msePerTrial(:,s,2)-r.msePerTrial(:,s,1);
    bootstrap = sort(mean(difference(randi(trials,trials,1000)),1));
    r.pairedMseDifferenceCI(s,:) = [mean(difference),bootstrap(25),bootstrap(975)];
end
end
