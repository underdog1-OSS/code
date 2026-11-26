%% 原 WD-RANSAC 与“一致集反馈抽样＋最终局部集WLS＋J择优”，打开后点击运行。
% 两种方法各评价3000个随机候选，前1000对点共用；改进方法末尾额外做一次WLS。
% 权重只在初始化最佳模型确定后和最佳J刷新时累加，不参与局部随机抽样。
% 末尾WLS只使用最终局部集中的点，保留这些点各自的累计权重。
% 有效WLS额外计算一次J；比本方法搜索结束时的最佳J严格更高才接受，平分保留原模型。
% 每个图窗采用2行2列，最多4个子图；超出时自动新建图窗。
% 显示原有对比、初始化及重估前后的汇总和逐次误差、后2000次最佳J刷新次数。
% 完整结果留在工作区result，不自动保存；真值不参与任何输出选择。
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(scriptDir),'modules'));

%% 1. 保留已确定的仿真设置。
cfg = wdr_config('paper');       % 每个SNR 1000次。
cfg.snrDb = -15:5:5;             % -15、-10、-5、0、5 dB。
cfg.mList = 3000;
cfg.initialBudget = 1000;        % 包含在3000次内，后续只有2000次。
cfg.consensusTauBins = 3;        % 固定一致集门限tau=3*实际DeltaOmega。
% 当前DeltaOmega=pi，tau=3*pi；单位沿用原代码的角频率/时间单位。
% 一致集门限保持固定；累计权重的三个门限按论文规则动态调整。
% 搜索中最佳模型始终按J判优；权重只用于末尾最终局部集的多点重估。
% 全局与局部集合中的同一点共享累计权重；退出局部集时不清除历史分数。
% 基线仍直接调用原wdr_estimate，使用原高斯抽样和分块评分。

%% 2. 一次预热，用于正式计时；不运行验证模块或参数试跑。
warmCfg = cfg; warmCfg.trials=1; warmCfg.snrDb=-5;
warmCfg.seed = cfg.seed+10000000; warmCfg.showProgress=false;
warmup = wdr_consensus_compare(warmCfg);
warmupSeconds = warmup.meanTotalSeconds;
clear warmup;
fprintf('正式统计：%d个SNR，每点%d次；M=%d，B=%d。\n', ...
    numel(cfg.snrDb),cfg.trials,cfg.mList,cfg.initialBudget);

%% 3. 正式统计：所有试验参加直接误差平均，真值只用于离线评价。
result = wdr_consensus_compare(cfg);
result.warmupSeconds = warmupSeconds;
% result.searchBestScore：随机搜索中实际计算得到的最高J。
% result.weightedRefitScore：原始WLS候选的J；result.score：择优后实际输出的J。
% 以下四项的第1维是试验编号，第2维是SNR，保存在线性域，画图时才转dB：
% result.initialError：前1000次初始化结束时，最高J模型的误差。
% result.preRefitError：本方法完成3000次搜索、尚未加权时的误差。
% result.postRefitError：原始WLS候选的误差，即使它被J择优拒绝也保留。
% result.error(:,:,2)：J择优后的最终误差；上述误差仅用于离线诊断。
% result.initialA、initialB：每次试验的初始模型参数；initialMseDb：初始化汇总MSE。
% result.weightedRefitApplied：有效WLS已完成；weightedRefitAccepted：最终选择了WLS。
% result.additionalScoreCount(:,:,2)：有效WLS时为1，其余为0，随机候选始终3000个。
% result.finalWeights：每次试验全部256个点的最终累计权重。
% result.refitWeights：末尾WLS实际使用的权重，最终局部集之外为零。
% result.refitPositiveWeightCount：末尾WLS实际参与拟合的正权重点数。
% result.examples：每个SNR的一次完整候选、计分门限与加分历史。

%% 4. 每个图窗最多4个子图，当前9个子图自动分成3个图窗（4+4+1）。
% 逐次误差横轴是Monte Carlo试验编号，不是单次搜索中的3000个候选编号。
% 图窗1为4个汇总图；后续图窗按SNR排列逐次误差，每窗最多4个。
% 手动保存图片建议按图窗编号命名：WD_RANSAC_Diagnostics_1.svg等。
nSnr = numel(result.snrDb);
nPlot = nSnr+4;
nFigure = ceil(nPlot/4);
figures = gobjects(nFigure,1);
layouts = cell(nFigure,1);
plotAxes = gobjects(nPlot,1);
for page = 1:nFigure
    if page==1
        pageName = '汇总对比';
    else
        pageName = '各SNR逐次误差';
    end
    figures(page) = figure('Name',sprintf('WD-RANSAC诊断 %d/%d：%s', ...
        page,nFigure,pageName),'Color','w', ...
        'Position',[40+25*(page-1) 40+25*(page-1) 1350 800]);
    layouts{page} = tiledlayout(figures(page),2,2, ...
        'TileSpacing','compact','Padding','compact');
    firstPlot = (page-1)*4+1;
    lastPlot = min(page*4,nPlot);
    for panel = firstPlot:lastPlot
        plotAxes(panel) = nexttile(layouts{page});
    end
    title(layouts{page},sprintf( ...
        '第%d/%d窗；M=%d；B=%d；tau=%.6g；最终局部集WLS＋J择优；每个SNR %d次', ...
        page,nFigure,cfg.mList,cfg.initialBudget,result.tau,cfg.trials));
end
initialColor = [0.6 0.15 0.75];
beforeColor = [0 0.35 0.85];
weightedColor = [0.9 0.2 0.1];
selectedColor = [0.1 0.55 0.15];
diagnosticNames = {sprintf('初始化结束（%d次）',cfg.initialBudget), ...
    sprintf('搜索结束（%d次，未加权）',cfg.mList), ...
    '加权后（未择优）','J择优输出'};

% 原有两方法对比，改进方法曲线使用J择优后的实际输出。
ax = plotAxes(1);
plot(ax,result.snrDb,result.mseDb(:,1),'b-o','LineWidth',1.6); hold(ax,'on');
plot(ax,result.snrDb,result.mseDb(:,2),'r-s','LineWidth',1.6);
grid(ax,'on'); xlabel(ax,'SNR / dB（原论文口径）');
ylabel(ax,'最终 IF MSE / dB [(rad/时间单位)^2]');
title(ax,'两方法最终输出的MSE');
legend(ax,result.methodNames,'Location','best','FontSize',8);
ax = plotAxes(2);
plot(ax,result.snrDb,result.meanTotalSeconds(:,1),'b-o','LineWidth',1.6); hold(ax,'on');
plot(ax,result.snrDb,result.meanTotalSeconds(:,2),'r-s','LineWidth',1.6);
grid(ax,'on'); xlabel(ax,'SNR / dB（原论文口径）');
ylabel(ax,'单次平均总运行时间 / 秒');
title(ax,'完整耗时：含局部集WLS及J评分和择优');
legend(ax,result.methodNames,'Location','best','FontSize',8);

% 四个阶段都属于改进算法；先平均线性误差再转dB，与左上主图的口径一致。
ax = plotAxes(3);
plot(ax,result.snrDb,result.initialMseDb,'-.^','Color',initialColor,'LineWidth',1.3);
hold(ax,'on');
plot(ax,result.snrDb,result.preRefitMseDb,'-o','Color',beforeColor,'LineWidth',1.3);
plot(ax,result.snrDb,result.postRefitMseDb,'-s','Color',weightedColor,'LineWidth',1.3);
plot(ax,result.snrDb,result.mseDb(:,2),'-d','Color',selectedColor,'LineWidth',1.3);
grid(ax,'on'); xlabel(ax,'SNR / dB'); ylabel(ax,'IF MSE / dB');
title(ax,'初始模型与后续搜索、WLS、最终输出的误差');
legend(ax,diagnosticNames,'Location','best','FontSize',8);

% 每个SNR单独画全部试验，不平滑、不截断或排除大误差试验。
trialIndex = (1:cfg.trials).';
for s = 1:nSnr
    ax = plotAxes(4+s);
    plot(ax,trialIndex,10*log10(result.initialError(:,s)), ...
        '-.','Color',initialColor,'LineWidth',0.6);
    hold(ax,'on');
    plot(ax,trialIndex,10*log10(result.preRefitError(:,s)), ...
        '-','Color',beforeColor,'LineWidth',0.6);
    plot(ax,trialIndex,10*log10(result.postRefitError(:,s)), ...
        '-','Color',weightedColor,'LineWidth',0.6);
    plot(ax,trialIndex,10*log10(result.error(:,s,2)), ...
        '.','Color',selectedColor,'MarkerSize',4);
    grid(ax,'on'); xlabel(ax,'Monte Carlo试验编号');
    ylabel(ax,'单次 IF MSE / dB');
    xlim(ax,[0.5 cfg.trials+0.5]);
    title(ax,sprintf('SNR %+g dB；接受WLS %d/%d次', ...
        result.snrDb(s),sum(result.weightedRefitAccepted(:,s)),cfg.trials));
    legend(ax,diagnosticNames,'Location','best','FontSize',8);
end

% 记录本方法初始化之后的J纪录刷新次数，不包括初始B次和最终WLS。
ax = plotAxes(4);
plot(ax,result.snrDb,mean(result.updateCount,1).','-o', ...
    'Color',[0.5 0.15 0.7],'LineWidth',1.5);
grid(ax,'on'); xlabel(ax,'SNR / dB'); ylabel(ax,'平均最佳J刷新次数');
title(ax,sprintf('后%d次局部搜索的J纪录更新',cfg.mList-cfg.initialBudget));
disp(result.summary);
fprintf('\n完成：%d个子图已分为%d个图窗，每窗最多4个；结果保留在工作区result。\n', ...
    nPlot,nFigure);
fprintf('改进算法各SNR的异常全局回退总次数：\n');
disp(table(result.snrDb,sum(result.fallbackCount,1).', ...
    'VariableNames',{'SNR_dB','FallbackCount'}));
fprintf('改进算法初始化、搜索结束、原始WLS与J择优后的误差及接受情况：\n');
disp(table(result.snrDb,result.initialMseDb,result.preRefitMseDb, ...
    result.postRefitMseDb,result.mseDb(:,2), ...
    mean(result.weightedRefitAccepted,1).',sum(~result.weightedRefitApplied,1).', ...
    'VariableNames',{'SNR_dB','InitialMse_dB','BeforeWlsMse_dB','RawWlsMse_dB','SelectedMse_dB', ...
    'WlsAcceptedFraction','WlsNotAppliedCount'}));
fprintf('改进算法平均权重更新次数、有效WLS比例及额外J评分次数：\n');
disp(table(result.snrDb,mean(result.weightUpdateCount,1).', ...
    mean(result.weightedRefitApplied,1).',mean(result.additionalScoreCount(:,:,2),1).', ...
    'VariableNames',{'SNR_dB','MeanWeightUpdates','WeightedRefitFraction','MeanExtraJCount'}));
