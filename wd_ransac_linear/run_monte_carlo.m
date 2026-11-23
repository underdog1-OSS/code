%% 论文 Example 1：多次随机实验，绘制瞬时频率 MSE 随 SNR 变化的曲线
% 在编辑器中点击“运行”；第一次可把下面的 preset 改为 'quick'。
% 'paper'：每个 SNR 1000 次；'quick'：每个 SNR 50 次，用于检查流程。
% 本代码绘制 WD 峰值法和各个 M 下的原始 WD-RANSAC 基线。
% 原文没有交代的参数分布、高斯抽样参数等，已在 wdr_config.m 中注明。
% 运行中可在命令行窗口看到进度，无需在那里输入任何命令。
% 计算结果保留在工作区 result 变量中，不自动生成数据、图片或日志文件。
% Open and click Run. The full paper-sized experiment can take many minutes.
% To learn the code first, change preset to 'quick' (50 trials at each SNR).
% This implements WD maxima and WD-RANSAC, including M=10 from Fig. 1.
% See wdr_config.m for the explicitly documented choices where the paper
% omits details. Do not interpret this as the author's original source code.
% No MAT, CSV, figure, log, or other output files are automatically written.
scriptDir = fileparts(mfilename('fullpath')); % 本文件所在的 wd_ransac_linear 文件夹
moduleDir = fullfile(fileparts(scriptDir), 'modules'); % 旁边的 modules 文件夹
addpath(moduleDir);
preset = 'paper';            % 'paper': 1000 trials/SNR; 'quick': 50 trials/SNR
cfg = wdr_config(preset);

%% 1. Monte Carlo 实验：每次重新生成参数、噪声和随机点对
result = wdr_monte_carlo(cfg);

%% 2. 先对平方误差求平均，再转换成 dB；不能先取 dB 再平均
figure('Name', 'WD-RANSAC Example 1: Monte Carlo', 'Color', 'w');
plot(result.snrDb, result.wdMseDb, 'k-', 'LineWidth', 1.8); hold on;
markers = {'o','*','+','d','s'};
labels = cell(1, numel(cfg.mList)+1);
labels{1} = 'WD maxima (direct error)';
for k = 1:numel(cfg.mList)
    plot(result.snrDb, result.ransacMseDb(:,k), ...
        ['-' markers{1+mod(k-1,numel(markers))}], 'LineWidth', 1.2);
    labels{k+1} = sprintf('WD-RANSAC M = %d', cfg.mList(k));
end
grid on; xlabel('SNR (dB, paper convention)');
ylabel('IF MSE (dB, angular frequency)');
title(sprintf('Example 1 implementation, %d trials/SNR', cfg.trials));
legend(labels, 'Location', 'northeast');
xlim([min(cfg.snrDb) max(cfg.snrDb)]);
fprintf('\nDone. Per-trial errors and settings remain in workspace variable result.\n');
fprintf('Curves are WD + WD-RANSAC only; see wdr_config for assumptions.\n');

