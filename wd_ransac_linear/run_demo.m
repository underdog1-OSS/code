%% WD-RANSAC 线性调频实验：先运行这个文件，看一组信号的估计结果
% 在 MATLAB 编辑器中打开本文件，直接点击“运行”。
% 如弹出文件不在当前文件夹的提示，选择“更改文件夹”。
% 本入口自动加载 modules 文件夹，不需要手动输入路径命令。
% 这是模块化代码：配置 -> 信号生成 -> WD -> 峰值 -> 随机抽点与相干评分。
% 图像只显示，计算结果保留在工作区，不自动保存任何输出文件。
% 固定 a,b 便于观察过程；论文中的随机参数实验请运行 run_monte_carlo.m。
% Open this file and click Run. No toolbox or command-window setup is needed.
% Modules: config -> signal -> WD -> ridge -> random pairs -> coherent score.`
% Source: I. Djurovic, IEEE SPL 23(5), 2016, DOI 10.1109/LSP.2016.2551732.
% This demo uses fixed a,b for inspection; run_monte_carlo uses random a,b.
% Only figures/workspace variables are produced; no files are saved.
scriptDir = fileparts(mfilename('fullpath')); % 本文件所在的 wd_ransac_linear 文件夹
moduleDir = fullfile(fileparts(scriptDir), 'modules'); % 旁边的 modules 文件夹
addpath(moduleDir);                             % 相当于手动添加路径
cfg = wdr_config('paper');
rng(cfg.seed, 'twister');

%% 1. 单次演示参数：可以先只改信噪比 snrDb，观察估计结果
snrDb = -5;                  % PAPER SNR = 10*log10(A^2/sigma^2)
aTrue = 220;                 % Angular-frequency slope, rad/(time unit)^2
bTrue = -70;                 % Angular frequency at t=0, rad/(time unit)
mList = [10 100 300 900 3000];

%% 2. 调用各模块：512 点信号，中间 256 点的瞬时频率估计
signal = wdr_generate_lfm(cfg, snrDb, [aTrue bTrue]);
[wd, omegaGrid, tWd] = wdr_distribution(signal.x, cfg);
ridge = wdr_ridge(wd, omegaGrid);
estimate = wdr_estimate(signal.x, cfg.t, tWd, ridge, mList, cfg);
omegaTrue = signal.a*tWd + signal.b;
% 直接误差：用于检查收窄真实频率范围后边界突变是否仍存在。
wdMse = mean((ridge - omegaTrue).^2);
% 周期误差的对照版本；启用时需注释上面一行。
% wdPeriod = pi/cfg.dt;
% wdDifference = ridge - omegaTrue;
% wdDifference = mod(wdDifference+wdPeriod/2,wdPeriod)-wdPeriod/2;
% wdMse = mean(wdDifference.^2);
% WDRANSAC 使用直接频率误差。
ransacMse = mean((estimate.omega - omegaTrue).^2, 1);
% 周期误差的对照版本，启用时需注释上面一行。
% ransacPeriod = pi/cfg.dt;
% ransacDifference = estimate.omega - omegaTrue;
% ransacDifference = mod(ransacDifference+ransacPeriod/2,ransacPeriod)-ransacPeriod/2;
% ransacMse = mean(ransacDifference.^2, 1);

%% 3. 显示 WD 和频率轨迹：这里用角频率，不是 Hz
figure('Name', 'WD-RANSAC linear FM: single trial', 'Color', 'w');
subplot(1,2,1);
imagesc(tWd, omegaGrid, wd); axis xy; colorbar; hold on;
plot(tWd, omegaTrue, 'w--', 'LineWidth', 1.5);
xlabel('t'); ylabel('\omega (rad / time unit)');
title(sprintf('WD, paper SNR = %.1f dB', snrDb));

subplot(1,2,2);
plot(tWd, ridge, '.', 'Color', [0.65 0.65 0.65]); hold on;
plot(tWd, omegaTrue, 'k-', 'LineWidth', 2);
plot(tWd, estimate.omega(:,end), 'r--', 'LineWidth', 1.8);
grid on; xlabel('t'); ylabel('\omega (rad / time unit)');
title(sprintf('WD-RANSAC, M = %d', mList(end)));
legend('WD maxima', 'True IF', 'WD-RANSAC', 'Location', 'best');

fprintf('\nTrue a = %.6f, b = %.6f\n', signal.a, signal.b);
fprintf('WD maxima direct IF MSE = %.6g, %.3f dB\n', wdMse, 10*log10(wdMse));
summary = table(mList(:), estimate.a(:), estimate.b(:), ...
    estimate.score(:), ransacMse(:), 10*log10(ransacMse(:)), ...
    'VariableNames', {'M','a_est','b_est','CoherentScore','Direct_IF_MSE','Direct_IF_MSE_dB'});
disp(summary);
fprintf('Paper SNR %.1f dB equals total-complex-noise SNR %.4f dB.\n', ...
    snrDb, snrDb - 10*log10(2));

