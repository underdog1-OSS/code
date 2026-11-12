%% WD-RANSAC + 质量排序和渐进抽样：单次对照
% 在 MATLAB 编辑器中打开本文件，点击“运行”；不自动保存图片或数据。
% 原版高斯抽样与渐进抽样使用同一条信号、同一个 WD 和相同候选数量。
% 质量分数是卡片中的待验证方案，不是 PROSAC 原文的 WD 质量公式。
% 不使用 PROSAC 停止规则，也不改变原两点成线与相干评分。
scriptDir = fileparts(mfilename('fullpath'));
moduleDir = fullfile(fileparts(scriptDir),'modules');
addpath(moduleDir);
cfg = wdr_config('paper');
snrDb = -5;
aTrue = 220;
bTrue = -70;
mList = cfg.mList;
growthBudget = max(mList); % PROSAC 的 T_N；所有 M 共用同一前缀扩展日程

rng(cfg.seed,'twister');
signal = wdr_generate_lfm(cfg,snrDb,[aTrue bTrue]);
timer = tic;
[wd,omega,tWd] = wdr_distribution(signal.x,cfg);
wdSeconds = toc(timer);
ridge = wdr_ridge(wd,omega);
% 原始 WD-RANSAC：高斯抽样。函数内部仍采用原来的抽样和相干评分。
timer = tic;
estimateGaussian = wdr_estimate(signal.x,cfg.t,tWd,ridge,mList,cfg);
gaussianSeconds = toc(timer);
% 排序和渐进抽样使用独立、可复现的随机序列。
rng(mod(cfg.seed+200000000,2^32),'twister');
[estimateProsac,quality,schedule,timing] = wdr_prosac_estimate( ...
    signal.x,cfg.t,tWd,ridge,wd,omega,mList,cfg,growthBudget);

% 真实频率仅用于仿真评价，没有用于 quality 或抽样。
omegaTrue = signal.a*tWd+signal.b;
gaussianMse = mean((estimateGaussian.omega-omegaTrue).^2,1);
prosacMse = mean((estimateProsac.omega-omegaTrue).^2,1);
tolerance = 2*median(diff(omega)); % 仿真标签容差：两个频率格，待验证
correct = abs(ridge-omegaTrue)<=tolerance;
n10 = max(1,ceil(0.10*numel(ridge)));
n20 = max(1,ceil(0.20*numel(ridge)));
fprintf('\nPoint accuracy: all %.1f%%, top 10%% %.1f%%, top 20%% %.1f%%\n', ...
    100*mean(correct),100*mean(correct(quality.order(1:n10))), ...
    100*mean(correct(quality.order(1:n20))));
fprintf('Diagnostic tolerance = %.4g rad/time unit (2 frequency bins).\n',tolerance);
fprintf('WD %.3f s | Gaussian sampling+scoring %.3f s | quality %.3f s | PROSAC sampling %.3f s | scoring %.3f s\n', ...
    wdSeconds,gaussianSeconds,timing.qualitySeconds,timing.samplingSeconds,timing.scoringSeconds);
summary = table(mList(:),gaussianMse(:),prosacMse(:), ...
    estimateGaussian.score(:),estimateProsac.score(:), ...
    'VariableNames',{'M','Gaussian_MSE','PROSAC_MSE','Gaussian_J','PROSAC_J'});
disp(summary);

figure('Name','WD-RANSAC and PROSAC: single trial','Color','w');
subplot(1,3,1);
plot(tWd,ridge,'.','Color',[0.7 0.7 0.7]); hold on;
plot(tWd,omegaTrue,'k-','LineWidth',1.5);
plot(tWd,estimateGaussian.omega(:,end),'b--','LineWidth',1.3);
plot(tWd,estimateProsac.omega(:,end),'r--','LineWidth',1.3);
grid on; xlabel('t'); ylabel('\omega (rad / time unit)');
title(sprintf('IF estimates, M = %d',mList(end)));
legend('WD maxima','True IF','Gaussian','PROSAC','Location','best');
subplot(1,3,2);
plot(tWd,quality.q,'.-'); grid on; xlabel('t'); ylabel('Quality q');
title('Point quality for ranking');
subplot(1,3,3);
plot(1:numel(schedule.prefixSize),schedule.prefixSize); grid on;
xlabel('Candidate number'); ylabel('Ranked prefix size');
title(sprintf('PROSAC growth, T_N = %d',growthBudget));
