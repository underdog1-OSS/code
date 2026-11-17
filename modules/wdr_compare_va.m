function result = wdr_compare_va(cfg,lambda,epsilonP,growthBudget,slopeThreshold)
% 对相同信号比较原最大峰高斯抽样、VA 点均匀抽样和 VA 点 PROSAC。
% 三者使用相同 M 和原相干评分；真值仅用于事后评价。
mList=cfg.mList(:).'; nM=numel(mList); nSnr=numel(cfg.snrDb); nt=cfg.trials;
if nargin<5, slopeThreshold=200; end
validateattributes(slopeThreshold,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(mList,{'numeric'},{'integer','positive','finite','nonempty'});
validateattributes(nt,{'numeric'},{'integer','positive','scalar'});
result.cfg=cfg; result.lambda=lambda; result.eta=0; result.epsilonP=epsilonP;
result.growthBudget=growthBudget; result.mList=mList; result.snrDb=cfg.snrDb(:);
result.slopeThreshold=slopeThreshold;
names={'baseline','vaUniform','vaProsac'};
for k=1:3
    result.([names{k} 'Error'])=zeros(nt,nSnr,nM);
    result.([names{k} 'SamplingSeconds'])=zeros(nt,nSnr);
    result.([names{k} 'ScoringSeconds'])=zeros(nt,nSnr);
end
result.wdSeconds=zeros(nt,nSnr); result.ridgeSeconds=zeros(nt,nSnr);
result.vaSeconds=zeros(nt,nSnr); result.clippedFraction=zeros(nt,nSnr);
result.nonpositiveFraction=zeros(nt,nSnr); result.tieColumns=zeros(nt,nSnr);
result.qualityTieFraction=zeros(nt,nSnr);
result.correctAll=zeros(nt,nSnr); result.correctTop10=zeros(nt,nSnr);
result.correctTop20=zeros(nt,nSnr); result.trialSeed=zeros(nt,nSnr);
result.trueSlope=zeros(nt,nSnr); result.selectedDifferent=zeros(nt,nSnr);
result.correctTolerance=2*pi/(cfg.wdLength*cfg.dt);
start=tic;
for s=1:nSnr
    for trial=1:nt
        seed=cfg.seed+(s-1)*nt+trial; rng(seed,'twister');
        signal=wdr_generate_lfm(cfg,cfg.snrDb(s));
        timer=tic; [wd,omega,tWd]=wdr_distribution(signal.x,cfg); result.wdSeconds(trial,s)=toc(timer);
        timer=tic; ridge=wdr_ridge(wd,omega); result.ridgeSeconds(trial,s)=toc(timer);
        timer=tic; pairs=wdr_sample_pairs(numel(ridge),max(mList),cfg); result.baselineSamplingSeconds(trial,s)=toc(timer);
        timer=tic; baseline=wdr_estimate(signal.x,cfg.t,tWd,ridge,mList,cfg,pairs); result.baselineScoringSeconds(trial,s)=toc(timer);

        timer=tic; va=wdr_va_select(wd,omega,tWd,lambda,epsilonP); result.vaSeconds(trial,s)=toc(timer);
        result.nonpositiveFraction(trial,s)=va.nonpositiveFraction;
        result.clippedFraction(trial,s)=va.clippedFraction; result.tieColumns(trial,s)=va.tieColumns;
        result.qualityTieFraction(trial,s)=va.qualityTieFraction;
        result.selectedDifferent(trial,s)=mean(va.selectedOmega~=ridge);
        result.trueSlope(trial,s)=signal.a;
        rng(mod(seed+100000000,2^32),'twister');
        timer=tic; pairs=zeros(2,max(mList));
        for k=1:max(mList), pairs(:,k)=sort(randperm(numel(ridge),2)); end
        result.vaUniformSamplingSeconds(trial,s)=toc(timer);
        timer=tic; uniform=wdr_estimate(signal.x,cfg.t,tWd,va.selectedOmega,mList,cfg,pairs); result.vaUniformScoringSeconds(trial,s)=toc(timer);
        rng(mod(seed+200000000,2^32),'twister');
        timer=tic; pairs=wdr_prosac_pairs(va.rankOrder,max(mList),growthBudget); result.vaProsacSamplingSeconds(trial,s)=toc(timer);
        timer=tic; prosac=wdr_estimate(signal.x,cfg.t,tWd,va.selectedOmega,mList,cfg,pairs); result.vaProsacScoringSeconds(trial,s)=toc(timer);

        result.baselineError(trial,s,:)=reshape(mean((baseline.omega-signal.omegaTrue).^2,1),1,1,nM);
        result.vaUniformError(trial,s,:)=reshape(mean((uniform.omega-signal.omegaTrue).^2,1),1,1,nM);
        result.vaProsacError(trial,s,:)=reshape(mean((prosac.omega-signal.omegaTrue).^2,1),1,1,nM);
        correct=abs(va.selectedOmega-signal.omegaTrue)<=result.correctTolerance;
        n10=ceil(0.1*numel(ridge)); n20=ceil(0.2*numel(ridge));
        result.correctAll(trial,s)=mean(correct);
        result.correctTop10(trial,s)=mean(correct(va.rankOrder(1:n10)));
        result.correctTop20(trial,s)=mean(correct(va.rankOrder(1:n20)));
        result.trialSeed(trial,s)=seed;
        if cfg.showProgress && (mod(trial,10)==0 || trial==nt)
            fprintf('SNR %+g dB | %d/%d 次 | 用时 %.1f 秒\n',cfg.snrDb(s),trial,nt,toc(start));
        end
    end
end
for k=1:3
    name=names{k}; mse=reshape(mean(result.([name 'Error']),1),nSnr,nM);
    result.([name 'Mse'])=mse; result.([name 'MseDb'])=10*log10(mse);
end
% 总时间包括 WD、选点/排序、抽样、评分；不含真值诊断和绘图。
result.baselineTotalSeconds=result.wdSeconds+result.ridgeSeconds+result.baselineSamplingSeconds+result.baselineScoringSeconds;
result.vaUniformTotalSeconds=result.wdSeconds+result.vaSeconds+result.vaUniformSamplingSeconds+result.vaUniformScoringSeconds;
result.vaProsacTotalSeconds=result.wdSeconds+result.vaSeconds+result.vaProsacSamplingSeconds+result.vaProsacScoringSeconds;
% 按预先设定的斜率幅度分组，只用于事后诊断，不参与选点和排序。
result.lowSlopeCount=zeros(nSnr,1); result.highSlopeCount=zeros(nSnr,1);
result.lowSlopeTop10=nan(nSnr,1); result.highSlopeTop10=nan(nSnr,1);
for s=1:nSnr
    low=abs(result.trueSlope(:,s))<=slopeThreshold; high=~low;
    result.lowSlopeCount(s)=sum(low); result.highSlopeCount(s)=sum(high);
    if any(low), result.lowSlopeTop10(s)=mean(result.correctTop10(low,s)); end
    if any(high), result.highSlopeTop10(s)=mean(result.correctTop10(high,s)); end
end
result.elapsedSeconds=toc(start);
end
