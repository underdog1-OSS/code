%% 高斯基线、均匀抽样与 PROSAC 抽样：固定 M 的蒙特卡洛对照
% 在 MATLAB 编辑器中打开本文件，点击“运行”。首次建议使用 quick。
% 仅显示图像和工作区变量 result，不自动生成图片、数据或日志文件。
% 三种方法共用每次试验的信号、WD、峰点及 M；均使用原相干评分。
% 排序质量是待验证的工程假设；不使用 PROSAC 的自适应停止规则。
scriptDir = fileparts(mfilename('fullpath'));
moduleDir = fullfile(fileparts(scriptDir),'modules');
addpath(moduleDir);
preset = 'quick';                 % quick: 每个 SNR 50 次；paper: 1000 次
cfg = wdr_config(preset);
growthBudget = max(cfg.mList);    % T_N；所有 M 共用同一渐进日程
result = wdr_compare_prosac(cfg,growthBudget);

% 图只显示最大的 M；完整各 M 结果保存在 result 中。
figure('Name','WD-RANSAC sampling comparison','Color','w');
subplot(1,2,1);
plot(result.snrDb,result.gaussianMseDb(:,end),'b-o','LineWidth',1.2); hold on;
plot(result.snrDb,result.uniformMseDb(:,end),'k-s','LineWidth',1.2);
plot(result.snrDb,result.prosacMseDb(:,end),'r-^','LineWidth',1.2);
grid on; xlabel('SNR (dB, paper convention)'); ylabel('Direct IF MSE (dB)');
title(sprintf('M = %d, %d trials/SNR',cfg.mList(end),cfg.trials));
legend('Original Gaussian','Uniform','Quality + PROSAC','Location','best');
subplot(1,2,2);
plot(result.snrDb,mean(result.correctAll,1),'k-o'); hold on;
plot(result.snrDb,mean(result.correctTop10,1),'r-^');
plot(result.snrDb,mean(result.correctTop20,1),'b-s');
grid on; xlabel('SNR (dB, paper convention)'); ylabel('Fraction close to true IF');
title('Quality ranking diagnostic (simulation only)');
legend('All peaks','Top 10%','Top 20%','Location','best');

fprintf('\nQuality diagnostic tolerance: %.4g rad/time unit (2 frequency bins).\n', ...
    result.correctTolerance);
fprintf('Average seconds/trial: WD %.4f | quality %.4f | Gaussian sampling+scoring %.4f | uniform sampling+scoring %.4f | PROSAC sampling+scoring %.4f\n', ...
    mean(result.wdSeconds(:)),mean(result.qualitySeconds(:)), ...
    mean(result.gaussianSamplingSeconds(:)+result.gaussianScoringSeconds(:)), ...
    mean(result.uniformSamplingSeconds(:)+result.uniformScoringSeconds(:)), ...
    mean(result.prosacSamplingSeconds(:)+result.prosacScoringSeconds(:)));
fprintf('PROSAC total beyond WD includes both quality and sampling+scoring.\n');
