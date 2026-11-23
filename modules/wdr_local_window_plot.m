function wdr_local_window_plot(result)
% 只显示最终频率轨迹 MSE 和完整算法总耗时两张图，不自动保存文件。
cfg = result.cfg;
snr = result.snrDb;
colors = [0.1 0.35 0.85;0.85 0.2 0.1;0.1 0.6 0.25];
markers = {'o','s','^'};
names = cfg.methodNames;

% 图一：各次试验的轨迹 MSE 取平均，纵轴采用对数刻度，单位为 Hz^2。
% 如需手动保存，建议文件名：local_window_MSE_comparison。
fig = figure('Name','局部窗口取峰：最终 MSE 对比', ...
    'Color','w','Visible','on','Position',[100 100 900 600]);
ax = axes('Parent',fig);
hold(ax,'on');
for m = 1:3
    mseHz2 = mean(result.trackMseHz2(:,:,m),1);
    semilogy(ax,snr,max(mseHz2,realmin),['-' markers{m}], ...
        'Color',colors(m,:),'LineWidth',1.6,'MarkerSize',7);
end
set(ax,'YScale','log');
grid(ax,'on');
xlabel(ax,'SNR / dB（原论文口径）');
ylabel(ax,'最终频率轨迹 MSE / Hz^2');
title(ax,sprintf('最终 MSE 对比：M=%d，每个 SNR %d 次',cfg.mList,cfg.trials));
legend(ax,names,'Location','best','Interpreter','none');

% 图二：每个 SNR 下完整算法单次运行时间的均值，不拆分各步骤。
% 总时间包括 WD、取峰、抽样和相干评分，不含试跑、绘图及诊断统计。
% 如需手动保存，建议文件名：local_window_total_runtime_comparison。
fig = figure('Name','局部窗口取峰：总运行时间对比', ...
    'Color','w','Visible','on','Position',[150 150 900 600]);
ax = axes('Parent',fig);
hold(ax,'on');
for m = 1:3
    totalMilliseconds = 1000*mean(result.totalSeconds(:,:,m),1);
    plot(ax,snr,totalMilliseconds,['-' markers{m}], ...
        'Color',colors(m,:),'LineWidth',1.6,'MarkerSize',7);
end
grid(ax,'on');
xlabel(ax,'SNR / dB（原论文口径）');
ylabel(ax,'单次平均总运行时间 / 毫秒');
title(ax,sprintf('总运行时间对比：M=%d，每个 SNR %d 次',cfg.mList,cfg.trials));
legend(ax,names,'Location','best','Interpreter','none');
end
