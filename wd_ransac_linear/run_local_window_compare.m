%% 小规模三组对照：打开本脚本后点击“运行”即可完整复现。
% 局部窗口参数在 wdr_local_window_config.m，三组共用原 WD-RANSAC 后端。
% 运行后只显示最终 MSE 和总运行时间两张图；数据留在工作区，不自动保存文件。
scriptDir = fileparts(mfilename('fullpath'));
codeDir = fileparts(scriptDir);
addpath(scriptDir);
addpath(fullfile(codeDir,'modules'));
cfg = wdr_local_window_config;

%% 1. 加噪实验前的无噪声、边界和递推检查，失败时立即停止。
verification = verify_local_window;

%% 2. 用极少量试验确认耗时；不使用试跑结果挑选参数或筛掉失败试验。
pilotCfg = cfg;
pilotCfg.trials = cfg.pilotTrials;
pilotCfg.snrDb = cfg.pilotSnrDb;
pilotCfg.showProgress = false;
pilot = wdr_local_window_compare(pilotCfg);
fprintf('\n计时试跑：每个 SNR %d 次，三组单次平均 %.4f / %.4f / %.4f 秒。\n', ...
    pilotCfg.trials,mean(pilot.totalSeconds(:,:,1),'all'), ...
    mean(pilot.totalSeconds(:,:,2),'all'),mean(pilot.totalSeconds(:,:,3),'all'));
fprintf('正式规模固定为 %d 个 SNR × %d 次 × 3 组，不自动扩大。\n',numel(cfg.snrDb),cfg.trials);

%% 3. 正式对照：所有预定试验全部参与统计。
result = wdr_local_window_compare(cfg);
result.verification = verification;
result.pilot = pilot;

%% 4. 显示最终 MSE 和总运行时间对照，结果表保留在工作区供查看。
wdr_local_window_plot(result);
fprintf('\n实验完成，图窗已显示，数据保留在工作区 result。\n');
disp(result.summary(:,{'SNR_dB','Method','TrackMSE_Hz2','Total_Seconds'}));
