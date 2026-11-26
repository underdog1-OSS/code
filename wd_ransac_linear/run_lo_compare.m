%% 师兄WD-LO-RANSAC复现对照：打开本脚本并点击运行。
% 只生成一个大图：左侧最终IF MSE，右侧平均单次总运行时间。
% 不自动保存文件，不新建result文件夹；数值留在工作区result。
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(scriptDir),'modules'));

%% 1. 仿真配置；修改这里即可改变SNR、试验次数和随机候选预算。
cfg = wdr_config('paper');
cfg.snrDb = -15:1:5;       % 师兄稿件范围，单位dB；本入口采用复噪声总功率SNR。
cfg.trials = 1000;         % 正式试验次数可修改；稿件为10000次。
cfg.mList = 3000;          % 两方法共享全部3000个全局随机两点候选。
cfg.loTauBins = 3;         % 一致点门限为3个实际频率格点，即当前3*pi rad/s。
% 固定使用wdr_config中的原始种子，不因切换算法而改变测试样本。
% 两方法共同按abs(WD)取峰；旧基线函数和旧入口保持原来的取峰与SNR口径。
% 本入口SNR若记为S，旧入口对应SNR为S+3.0103dB，比较旧图时需换算。
% 每当原候选J超过当前全局最优才触发一次WLS；不做累计权重或局部抽样。

%% 2. 无噪声轴检查及一次预热；不产生图像或写入文件。
checkCfg = cfg; checkCfg.mList = 8;
clean = wdr_generate_lfm(checkCfg,Inf,[0 23*pi]);
[checkWd,checkOmega,checkTime] = wdr_distribution(clean.x,checkCfg);
[checkRidge,checkAmplitude] = wdr_ridge(abs(checkWd),checkOmega);
checkCfg.loTau = cfg.loTauBins*mean(diff(checkOmega));
assert(isequal(size(checkWd),[256 256]) && ...
    max(abs(checkTime-cfg.tWd))<1e-12 && max(abs(checkRidge-23*pi))<1e-10, ...
    '无噪声WD的维度、频率轴或时间轴检查失败。');
rng(cfg.seed+10000000,'twister');
checkPairs = wdr_sample_pairs(numel(checkRidge),checkCfg.mList,checkCfg);
checkEstimate = wdr_lo_estimate(clean.x,cfg.t,checkTime,checkRidge, ...
    checkAmplitude,checkCfg,checkPairs);
assert(max(abs(checkEstimate.omega-clean.omegaTrue))<1e-9 && ...
    abs(checkEstimate.score-cfg.nSamples*cfg.amplitude)<1e-8, ...
    '无噪声估计或全信号相干评分检查失败。');
clear clean checkWd checkOmega checkTime checkRidge checkAmplitude checkPairs checkEstimate checkCfg;
warmCfg = cfg; warmCfg.trials = 1; warmCfg.snrDb = 0;
warmCfg.seed = cfg.seed+10000000; warmCfg.showProgress = false;
warmup = wdr_lo_compare(warmCfg);
clear warmup warmCfg;

%% 3. 正式统计；算法执行顺序交替，真值误差和绘图不计入算法耗时。
fprintf('师兄算法对照：%d个SNR，每点%d次，M=%d，tau=%g*pi rad/s。\n', ...
    numel(cfg.snrDb),cfg.trials,cfg.mList,cfg.loTauBins);
fprintf('两方法统一使用abs(WD)取峰、复噪声总功率SNR，保留全部试验。\n');
result = wdr_lo_compare(cfg);
disp(result.summary);

%% 4. 只生成一个figure和两个子图，图片由你手动保存。
% 建议文件名：WD_RANSAC_vs_WD_LO_RANSAC_M3000_1000trials.svg。
figure('Name','WD-RANSAC与师兄WD-LO-RANSAC对照','Color','w', ...
    'Position',[80 100 1250 470]);
layout = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile;
plot(result.snrDb,result.mseDb(:,1),'b-o','LineWidth',1.5); hold on;
plot(result.snrDb,result.mseDb(:,2),'r-s','LineWidth',1.5);
grid on; xlabel('SNR / dB（复噪声总功率）');
ylabel('最终 IF MSE / dB [(rad/s)^2]');
title('最终估计误差'); legend(result.methodNames,'Location','best');
nexttile;
plot(result.snrDb,result.meanTotalSeconds(:,1),'b-o','LineWidth',1.5); hold on;
plot(result.snrDb,result.meanTotalSeconds(:,2),'r-s','LineWidth',1.5);
grid on; xlabel('SNR / dB（复噪声总功率）');
ylabel('平均单次总运行时间 / 秒');
title('包含WD、取峰、抽样和全部评分');
legend(result.methodNames,'Location','best');
title(layout,sprintf('M=%d；tau=%.6g rad/s；每个SNR %d次', ...
    cfg.mList,result.tau,cfg.trials));
