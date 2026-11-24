%% 原 WD-RANSAC 与最终最优直线精化对照：打开本脚本，点击“运行”。
% 先选出原算法最高分直线，只精化这一条，再比较精化前后的相干得分。
% 精化得分严格提高才采用新直线；相同或降低时保留原直线。
% 只显示最终 MSE 和平均总运行时间两张图；数据留在 result，不保存文件。
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(scriptDir),'modules'));

%% 1. 本次参数集中放在这里。
cfg = wdr_config('paper');
cfg.snrDb = [-15 -10 -7 -5 0 5 10]; % SNR 单位为 dB，上限延长到 10 dB。
cfg.trials = 500;                % 先小规模比较，之后可以自行改为 1000。
cfg.mList = 3000;                % 两组共用相同的 3000 对抽样时刻。
cfg.seed = 20261124;
cfg.osheaWindowLength = 17;      % 奇数平均长度；这是实现选择，不是原文规定值。
% 模型为 exp(j*(a*t^2/2+b*t))，时间沿用基线的归一化时间单位。
% a 为 rad/时间单位^2，b 为 rad/时间单位；MSE 使用直接角频率误差。
% SNR=10*log10(A^2/sigma^2)，噪声实部、虚部分别有方差 sigma^2。
% 精化步骤沿用原文的去啁啾、平均、展开和重估，接入位置改为最终输出后。
% 这是本次要求的后处理版本；单个 M 的每次试验只调用一次精化。
% 本实现的移动平均仅使用完整居中窗口；17 点时用中间 496 点拟合相位。
% 粗估计离真值太远时，残余相位展开可能失败，精化不保证一定有效。

%% 2. 无噪声小检查随本脚本运行，不出图，也不计入两组算法耗时。
check_oshea_module(cfg);

%% 3. 相同信号、WD、峰点和抽样时刻下进行配对比较。
result = wdr_compare_oshea(cfg);

%% 4. 只显示最终频率轨迹 MSE。
% 若手动保存，建议文件名：WD_RANSAC_OShea_FinalOnly_MSE_vs_SNR.svg。
fig = figure('Name','WD-RANSAC 最优直线精化择优：最终 MSE', ...
    'Color','w','Position',[100 100 900 600]);
ax = axes('Parent',fig); hold(ax,'on');
plot(ax,result.snrDb,result.mseDb(:,1),'b-o','LineWidth',1.6,'MarkerSize',7);
plot(ax,result.snrDb,result.mseDb(:,2),'r-s','LineWidth',1.6,'MarkerSize',7);
grid(ax,'on');
xlabel(ax,'SNR / dB（原 WD-RANSAC 论文口径）');
ylabel(ax,'频率轨迹 MSE / dB（角频率平方误差）');
title(ax,sprintf('M=%d，每个 SNR %d 次；移动平均 %d 点', ...
    cfg.mList,cfg.trials,cfg.osheaWindowLength));
legend(ax,result.methodNames,'Location','best','Interpreter','none');

%% 5. 只显示每个 SNR 下完整算法的单次平均总耗时。
% 若手动保存，建议文件名：WD_RANSAC_OShea_FinalOnly_Runtime_vs_SNR.svg。
fig = figure('Name','WD-RANSAC 最优直线精化择优：总运行时间', ...
    'Color','w','Position',[150 150 900 600]);
ax = axes('Parent',fig); hold(ax,'on');
plot(ax,result.snrDb,result.meanTotalMilliseconds(:,1),'b-o','LineWidth',1.6,'MarkerSize',7);
plot(ax,result.snrDb,result.meanTotalMilliseconds(:,2),'r-s','LineWidth',1.6,'MarkerSize',7);
grid(ax,'on');
xlabel(ax,'SNR / dB（原 WD-RANSAC 论文口径）');
ylabel(ax,'单次平均总运行时间 / 毫秒');
title(ax,sprintf('完整流程总耗时：M=%d，每个 SNR %d 次',cfg.mList,cfg.trials));
legend(ax,result.methodNames,'Location','best','Interpreter','none');
disp(result.summary);
fprintf('\n对照完成：两张图已显示，全部数据保留在工作区 result。\n');

function check_oshea_module(cfg)
% 用可解析的已知信号检查系数符号、二次相位因子和平均窗口的时间对齐。
rngState = rng;
cleaner = onCleanup(@() rng(rngState)); %#ok<NASGU>
rng(cfg.seed,'twister');
signal = wdr_generate_lfm(cfg,Inf,[220 -70]);
% 长度为 1 时不平滑，二次残余相位应能精确还原两个非恒定系数。
[a,b] = wdr_oshea_refine(signal.x,cfg.t,signal.a+3,signal.b-0.4,1);
assert(abs(a-signal.a)<1e-8 && abs(b-signal.b)<1e-8, ...
    'O''Shea 系数增量的符号或二次相位因子错误。');
% 只有残余频率偏差时，居中平均不能引入频率偏差或时间偏移。
[a,b,detail] = wdr_oshea_refine( ...
    signal.x,cfg.t,signal.a,signal.b+0.4,cfg.osheaWindowLength);
assert(abs(a-signal.a)<1e-8 && abs(b-signal.b)<1e-8, ...
    '移动平均的时间对齐或残余相位拟合错误。');
assert(numel(detail.fitIndices)==cfg.nSamples-cfg.osheaWindowLength+1);
% 相同的理想直线和点对输入时，两种评分都应达到完全相干幅值。
pairs = [1 40 100;256 210 180];
baseline = wdr_estimate(signal.x,cfg.t,cfg.tWd,signal.omegaTrue,3,cfg,pairs);
refined = wdr_oshea_estimate(signal.x,cfg.t,cfg.tWd,signal.omegaTrue,3,cfg,pairs);
assert(abs(baseline.score-cfg.nSamples*cfg.amplitude)<1e-8);
assert(abs(refined.score-baseline.score)<1e-8);
assert(max(abs(baseline.candidateA-refined.coarseCandidateA))<1e-10 && ...
    max(abs(baseline.candidateB-refined.coarseCandidateB))<1e-10);
assert(refined.refinementCount==1,'单个 M 只能精化一条最终最高分直线。');
% 含噪小例检查：精化来源必须是各预算下原算法的最高分直线。
noisy = wdr_generate_lfm(cfg,-7,[220 -70]);
[wd,omega,tWd] = wdr_distribution(noisy.x,cfg);
ridge = wdr_ridge(wd,omega);
before = wdr_estimate(noisy.x,cfg.t,tWd,ridge,[1 3],cfg,pairs);
after = wdr_oshea_estimate(noisy.x,cfg.t,tWd,ridge,[1 3],cfg,pairs);
assert(isequal(after.bestIteration,before.bestIteration) && ...
    isequal(after.bestPairs,before.bestPairs), ...
    '精化来源必须是原 WD-RANSAC 已经选出的最高分直线。');
assert(after.refinementCount==numel(before.mList), ...
    '每个 M 预算只能精化一条胜出直线。');
assert(all(after.score>=before.score-1e-9), ...
    '最终评分择优后，得分不能低于原最高分直线。');
% 独立重算最终输出模型的分数，避免分数与返回参数不对应。
finalPhase = (cfg.t.^2/2)*after.a+cfg.t*after.b;
manualScore = abs(noisy.x.'*exp(-1j*finalPhase));
assert(max(abs(manualScore-after.score))<1e-9, ...
    '最终输出直线的参数与保存的相干得分不一致。');
fprintf('检查通过：精化系数、窗口对齐、只精化原最高分直线及前后评分择优。\n');
end
