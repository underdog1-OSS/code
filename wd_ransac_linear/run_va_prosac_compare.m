%% VA 选点与 PROSAC 对照：在 MATLAB 编辑器中直接点击运行
% 不修改原 WD、旧质量分数和现有脚本，也不自动保存图片或数据。
% 建议图片名：WD_VA_PROSAC_comparison_lambda0p032_M3000_10trials.svg
scriptDir=fileparts(mfilename('fullpath')); addpath(fullfile(fileparts(scriptDir),'modules'));
cfgVA=wdr_config('quick');
cfgVA.trials=1000;                         % 初步检查；正式统计应增加试验次数
cfgVA.snrDb=[-9 -7 -5 -3 0];            % 初步检查的信噪比
cfgVA.mList=[10 100 300 900 3000];
lambdaVA=0.032;                          % 待验证的初始值，频率格编号差平方的系数
epsilonP=1e-12;                          % 全图归一化后的无量纲对数下限
slopeThreshold=200;                      % 仿真事后分组，角频率斜率单位／时间
growthBudget=max(cfgVA.mList);           % 固定日程，不使用概率停止规则
vaResult=wdr_compare_va(cfgVA,lambdaVA,epsilonP,growthBudget,slopeThreshold);

figure('Name','VA 选点与渐进抽样对照','Color','w','Position',[100 100 1300 780]);
subplot(2,2,1);
plot(vaResult.snrDb,vaResult.baselineMseDb(:,end),'b-o','LineWidth',1.2); hold on;
plot(vaResult.snrDb,vaResult.vaUniformMseDb(:,end),'k-s','LineWidth',1.2);
plot(vaResult.snrDb,vaResult.vaProsacMseDb(:,end),'r-^','LineWidth',1.2);
grid on; xlabel('信噪比 / dB'); ylabel('直接频率误差的 MSE / dB');
title(sprintf('M=%d；每个信噪比 %d 次',cfgVA.mList(end),cfgVA.trials));
legend('原最大峰＋高斯','VA 点＋均匀','VA 点＋PROSAC','Location','best');
subplot(2,2,2);
plot(vaResult.snrDb,mean(vaResult.correctAll,1),'k-o'); hold on;
plot(vaResult.snrDb,mean(vaResult.correctTop10,1),'r-^');
plot(vaResult.snrDb,mean(vaResult.correctTop20,1),'b-s');
grid on; xlabel('信噪比 / dB'); ylabel('所选点接近真实频率的比例');
title('仿真诊断：容差为两个频率格'); legend('全部点','前 10%','前 20%','Location','best');
subplot(2,2,3);
plot(vaResult.snrDb,mean(vaResult.baselineTotalSeconds,1),'b-o'); hold on;
plot(vaResult.snrDb,mean(vaResult.vaUniformTotalSeconds,1),'k-s');
plot(vaResult.snrDb,mean(vaResult.vaProsacTotalSeconds,1),'r-^');
grid on; xlabel('信噪比 / dB'); ylabel('每次试验平均时间 / 秒');
title('总时间：包含 WD、选点、抽样和评分'); legend('原最大峰＋高斯','VA 点＋均匀','VA 点＋PROSAC','Location','best');
subplot(2,2,4);
plot(vaResult.snrDb,mean(vaResult.selectedDifferent,1),'m-o');
grid on; xlabel('信噪比 / dB'); ylabel('VA 选点与最大峰不同的比例'); ylim([0 1]);
title(sprintf('lambda=%.3g，eta=0，epsilon=%.1e',lambdaVA,epsilonP));

fprintf('\n频率残差诊断容差：%.5g；lambda=%.5g，eta=0，epsilon=%.1e。\n', ...
    vaResult.correctTolerance,lambdaVA,epsilonP);
fprintf('每个信噪比的结果：M=%d，MSE 单位为 dB；时间为秒。\n',cfgVA.mList(end));
for s=1:numel(cfgVA.snrDb)
    fprintf('SNR=%+g：原基线 %.2f，VA均匀 %.2f，VA渐进 %.2f；VA计算 %.4f；总时间 %.4f / %.4f / %.4f\n', ...
        cfgVA.snrDb(s),vaResult.baselineMseDb(s,end),vaResult.vaUniformMseDb(s,end),vaResult.vaProsacMseDb(s,end), ...
        mean(vaResult.vaSeconds(:,s)),mean(vaResult.baselineTotalSeconds(:,s)), ...
        mean(vaResult.vaUniformTotalSeconds(:,s)),mean(vaResult.vaProsacTotalSeconds(:,s)));
    fprintf('  WD %.4f；原基线抽样 %.4f、评分 %.4f；VA均匀抽样 %.4f、评分 %.4f；VA渐进抽样 %.4f、评分 %.4f 秒。\n', ...
        mean(vaResult.wdSeconds(:,s)),mean(vaResult.baselineSamplingSeconds(:,s)),mean(vaResult.baselineScoringSeconds(:,s)), ...
        mean(vaResult.vaUniformSamplingSeconds(:,s)),mean(vaResult.vaUniformScoringSeconds(:,s)), ...
        mean(vaResult.vaProsacSamplingSeconds(:,s)),mean(vaResult.vaProsacScoringSeconds(:,s)));
end
fprintf('非正 WD 占比 %.2f%%；对数下限截断占比 %.2f%%；每次平均并列时刻 %.2f；排序分数重复占比 %.2f%%。\n', ...
    100*mean(vaResult.nonpositiveFraction(:)),100*mean(vaResult.clippedFraction(:)), ...
    mean(vaResult.tieColumns(:)),100*mean(vaResult.qualityTieFraction(:)));
fprintf('按斜率绝对值 %.1f 分组，仅用于仿真诊断；各组样本数少时不宜下结论。\n',slopeThreshold);
for s=1:numel(cfgVA.snrDb)
    fprintf('SNR=%+g：小斜率 %d 次，前10%%正确比例 %.3f；大斜率 %d 次，比例 %.3f。\n', ...
        cfgVA.snrDb(s),vaResult.lowSlopeCount(s),vaResult.lowSlopeTop10(s), ...
        vaResult.highSlopeCount(s),vaResult.highSlopeTop10(s));
end
