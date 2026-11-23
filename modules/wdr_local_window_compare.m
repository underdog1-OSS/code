function result = wdr_local_window_compare(cfg)
% 共用信号、WD、时间点对和相干评分，完整保留所有成功或失败试验。
nSnr = numel(cfg.snrDb); nTrials = cfg.trials;
nTime = numel(cfg.centerIndices); nMethods = 3;
assert(isscalar(cfg.mList), '本对照只使用一个共同的 M。');
result.cfg = cfg;
result.snrDb = cfg.snrDb(:);
result.timeSeconds = cfg.physicalTime(cfg.centerIndices);
fields = {'peakMseHz2','correctFraction','initialCorrect','windowTruthFraction', ...
    'maxErrorRun','errorRunCount','longErrorRun','f0ErrorHz','rateErrorHzPerSec', ...
    'trackMseHz2','f0HatHz','rateHatHzPerSec','score','bestIteration', ...
    'peakSeconds','fitScoreSeconds','totalSeconds','windowGivenCorrectPrevious', ...
    'errorAfterErrorFraction'};
for k = 1:numel(fields)
    result.(fields{k}) = nan(nTrials,nSnr,nMethods);
end
result.peakHz = zeros(nTime,nTrials,nSnr,nMethods);
result.windowLowerHz = nan(nTime,nTrials,nSnr,nMethods);
result.windowUpperHz = nan(nTime,nTrials,nSnr,nMethods);
result.trueF0Hz = zeros(nTrials,nSnr);
result.trueRateHzPerSec = zeros(nTrials,nSnr);
result.trialSeed = zeros(nTrials,nSnr);
result.wdSeconds = zeros(nTrials,nSnr);
result.pairSeconds = zeros(nTrials,nSnr);
result.methodOrder = zeros(nTrials,nSnr,nMethods);
clockStart = tic;
for s = 1:nSnr
    for trial = 1:nTrials
        % 固定间隔保证改变试验次数不会改变已有试验的随机种子。
        seed = cfg.seed+(s-1)*100000+trial;
        rng(seed,'twister');
        signal = wdr_local_window_signal(cfg,cfg.snrDb(s));
        timer = tic;
        [wd,omega,tWd] = wdr_distribution(signal.x,cfg);
        result.wdSeconds(trial,s) = toc(timer);
        timer = tic;
        pairs = wdr_sample_pairs(nTime,cfg.mList,cfg);
        result.pairSeconds(trial,s) = toc(timer);
        result.trialSeed(trial,s) = seed;
        result.trueF0Hz(trial,s) = signal.f0Hz;
        result.trueRateHzPerSec(trial,s) = signal.rateHzPerSec;
        % 平衡执行先后顺序，减小缓存、预热和后台负载的顺序偏差。
        order = circshift(1:3,[0 mod(trial+s-2,3)]);
        result.methodOrder(trial,s,:) = reshape(order,1,1,3);
        for m = order
            timer = tic;
            [ridge,detail] = wdr_local_window_ridge(wd,omega,cfg,char('A'+m-1),signal.omegaTrue(1));
            result.peakSeconds(trial,s,m) = toc(timer);
            timer = tic;
            estimate = wdr_estimate(signal.x,cfg.t,tWd,ridge,cfg.mList,cfg,pairs);
            result.fitScoreSeconds(trial,s,m) = toc(timer);
            result.totalSeconds(trial,s,m) = result.wdSeconds(trial,s)+ ...
                result.pairSeconds(trial,s)+result.peakSeconds(trial,s,m)+ ...
                result.fitScoreSeconds(trial,s,m);
            peakHz = ridge/(2*pi);
            correct = abs(peakHz-signal.trueHz)<=cfg.correctToleranceHz+1e-10;
            errors = ~correct;
            edges = diff([false;errors;false]);
            runLengths = find(edges==-1)-find(edges==1);
            result.maxErrorRun(trial,s,m) = max([0;runLengths]);
            result.errorRunCount(trial,s,m) = numel(runLengths);
            result.longErrorRun(trial,s,m) = any(runLengths>=cfg.longErrorRunLength);
            result.peakMseHz2(trial,s,m) = mean((peakHz-signal.trueHz).^2);
            result.correctFraction(trial,s,m) = mean(correct);
            result.initialCorrect(trial,s,m) = correct(1);
            if any(errors(1:end-1))
                result.errorAfterErrorFraction(trial,s,m) = ...
                    mean(errors(find(errors(1:end-1))+1));
            end
            if m>1
                contained = signal.omegaTrue>=detail.lower-1e-10 & ...
                    signal.omegaTrue<=detail.upper+1e-10;
                result.windowTruthFraction(trial,s,m) = mean(contained(2:end));
                previousCorrect = correct(1:end-1);
                if any(previousCorrect)
                    nextContained = contained(2:end);
                    result.windowGivenCorrectPrevious(trial,s,m) = mean(nextContained(previousCorrect));
                end
            end
            rateHat = estimate.a/(2*pi);
            f0Hat = (estimate.b+estimate.a*cfg.t(1))/(2*pi);
            result.rateHatHzPerSec(trial,s,m) = rateHat;
            result.f0HatHz(trial,s,m) = f0Hat;
            result.rateErrorHzPerSec(trial,s,m) = rateHat-signal.rateHzPerSec;
            result.f0ErrorHz(trial,s,m) = f0Hat-signal.f0Hz;
            result.trackMseHz2(trial,s,m) = mean(((estimate.omega-signal.omegaTrue)/(2*pi)).^2);
            result.score(trial,s,m) = estimate.score;
            result.bestIteration(trial,s,m) = estimate.bestIteration;
            result.peakHz(:,trial,s,m) = peakHz;
            result.windowLowerHz(:,trial,s,m) = detail.lower/(2*pi);
            result.windowUpperHz(:,trial,s,m) = detail.upper/(2*pi);
        end
        if cfg.showProgress && (mod(trial,20)==0 || trial==nTrials)
            fprintf('局部窗口对照：SNR %+g dB，%d/%d 次，累计 %.1f 秒\n', ...
                cfg.snrDb(s),trial,nTrials,toc(clockStart));
        end
    end
end
result.elapsedSeconds = toc(clockStart);
result.environment.matlabVersion = version;
result.environment.computer = computer;
result.environment.timeStamp = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
% 配对自助区间以“试验”为重抽样单位，不能把相邻 256 个时刻当独立样本。
result.paired = paired_summary(result);
result.summary = summary_table(result);
result.trials = trial_table(result);
end

function output = paired_summary(result)
cfg = result.cfg;
rngState = rng;
cleaner = onCleanup(@() rng(rngState)); %#ok<NASGU>
rng(cfg.bootstrapSeed,'twister');
metrics = {'correctFraction','trackMseHz2','maxErrorRun'};
contrasts = [2 1;3 1;2 3];
output = table;
for s = 1:numel(cfg.snrDb)
    for c = 1:size(contrasts,1)
        for k = 1:numel(metrics)
            x = result.(metrics{k})(:,s,contrasts(c,1))-result.(metrics{k})(:,s,contrasts(c,2));
            indexes = randi(numel(x),numel(x),cfg.bootstrapCount);
            distribution = sort(mean(x(indexes),1));
            lo = distribution(max(1,ceil(0.025*cfg.bootstrapCount)));
            hi = distribution(min(cfg.bootstrapCount,ceil(0.975*cfg.bootstrapCount)));
            contrast = string(char('A'+contrasts(c,1)-1))+'-'+string(char('A'+contrasts(c,2)-1));
            row = table(cfg.snrDb(s),contrast,string(metrics{k}),mean(x),lo,hi, ...
                'VariableNames',{'SNR_dB','Contrast','Metric','MeanDifference','CI95Low','CI95High'});
            output = [output;row]; %#ok<AGROW>
        end
    end
end
end

function output = summary_table(r)
output = table;
for s = 1:numel(r.snrDb)
    for m = 1:3
        row = table(r.snrDb(s),string(r.cfg.methodNames{m}),r.cfg.trials, ...
            mean(r.correctFraction(:,s,m)),mean(r.peakMseHz2(:,s,m)), ...
            mean(r.initialCorrect(:,s,m)),mean(r.windowTruthFraction(:,s,m),'omitnan'), ...
            mean(r.maxErrorRun(:,s,m)),max(r.maxErrorRun(:,s,m)), ...
            mean(r.longErrorRun(:,s,m)),mean(r.errorRunCount(:,s,m)), ...
            mean(r.windowGivenCorrectPrevious(:,s,m),'omitnan'), ...
            mean(r.errorAfterErrorFraction(:,s,m),'omitnan'), ...
            mean(r.trackMseHz2(:,s,m)),sqrt(mean(r.f0ErrorHz(:,s,m).^2)), ...
            sqrt(mean(r.rateErrorHzPerSec(:,s,m).^2)),mean(r.f0ErrorHz(:,s,m)), ...
            mean(r.rateErrorHzPerSec(:,s,m)),mean(r.wdSeconds(:,s)), ...
            mean(r.pairSeconds(:,s)),mean(r.peakSeconds(:,s,m)), ...
            mean(r.fitScoreSeconds(:,s,m)),mean(r.totalSeconds(:,s,m)), ...
            'VariableNames',{'SNR_dB','Method','Trials','PeakCorrectFraction','PeakMSE_Hz2', ...
            'InitialCorrectFraction','TruthInWindowFraction','MeanLongestErrorRun','WorstErrorRun', ...
            'LongRunTrialFraction','MeanErrorRunCount','TruthInWindowGivenCorrectPrevious', ...
            'ErrorAfterErrorFraction','TrackMSE_Hz2','F0_RMSE_Hz','Rate_RMSE_HzPerSec', ...
            'F0_Bias_Hz','Rate_Bias_HzPerSec','WD_Seconds','Pair_Seconds','Peak_Seconds', ...
            'FitScore_Seconds','Total_Seconds'});
        output = [output;row]; %#ok<AGROW>
    end
end
end

function output = trial_table(r)
output = table;
for s = 1:numel(r.snrDb)
    for m = 1:3
        n = r.cfg.trials;
        row = table(repmat(r.snrDb(s),n,1),(1:n)',r.trialSeed(:,s), ...
            repmat(string(r.cfg.methodNames{m}),n,1),r.trueF0Hz(:,s),r.trueRateHzPerSec(:,s), ...
            r.correctFraction(:,s,m),r.initialCorrect(:,s,m),r.windowTruthFraction(:,s,m), ...
            r.maxErrorRun(:,s,m),r.peakMseHz2(:,s,m),r.f0HatHz(:,s,m),r.rateHatHzPerSec(:,s,m), ...
            r.f0ErrorHz(:,s,m),r.rateErrorHzPerSec(:,s,m),r.trackMseHz2(:,s,m),r.totalSeconds(:,s,m), ...
            'VariableNames',{'SNR_dB','Trial','Seed','Method','TrueF0_Hz','TrueRate_HzPerSec', ...
            'PeakCorrectFraction','InitialCorrect','TruthInWindowFraction','LongestErrorRun', ...
            'PeakMSE_Hz2','F0Hat_Hz','RateHat_HzPerSec','F0Error_Hz','RateError_HzPerSec', ...
            'TrackMSE_Hz2','Total_Seconds'});
        output = [output;row]; %#ok<AGROW>
    end
end
end
