function cfg = wdr_local_window_config
% 局部取峰对照的独立配置；不修改原论文复现配置。
cfg = wdr_config('paper');
cfg.snrDb = [-15 -10 -7 -5 0 5];
cfg.trials = 100;
cfg.mList = 3000;
cfg.seed = 20261123;
cfg.showProgress = true;
% 本实验明确把原代码的时间单位解释为秒，保留原采样间隔和 WD。
cfg.fsHz = 1/cfg.dt;
cfg.physicalTime = cfg.t-cfg.t(1);
% 每次独立生成正频率上升调频信号；整个 512 点内也不接近频带边缘。
cfg.startFrequencyRangeHz = [16 28];
cfg.chirpRateRangeHzPerSec = [4 10];
cfg.windowFactors = [0.9 1.1];
cfg.correctToleranceBins = 1;
cfg.frequencyStepHz = 1/(2*cfg.wdLength*cfg.dt);
cfg.correctToleranceHz = cfg.correctToleranceBins*cfg.frequencyStepHz;
cfg.longErrorRunLength = 32;
cfg.methodNames = {'A_全频最大峰','B_实际初始化局部峰','C_理想初始化局部峰'};
cfg.bootstrapCount = 2000;
cfg.bootstrapSeed = 20261124;
cfg.pilotTrials = 3;
cfg.pilotSnrDb = [-10 0];
cfg.snrDefinition = '10log10(A^2/sigma^2)，噪声实部和虚部各有方差sigma^2；总复噪声功率口径比本值低3.0103dB';
cfg.errorDefinition = '所有组统一使用直接相减，不折回频率周期；频率为Hz，斜率为Hz/s';
cfg.initializationC = '首列使用真实频率最近的WD格点，之后只用上一列实际所选格点递推';
cfg.timingDefinition = '单方法时间=共用WD时间+该方法取峰时间+共用抽样时间+该方法拟合评分时间；不含信号生成、诊断和绘图';
end
