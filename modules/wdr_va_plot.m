function fig = wdr_va_plot(result)
% 只画一张大图的两个子图；在编辑器直接运行本函数可选择 MAT 文件重绘。
loadedFromFile = false;
if nargin<1
    [name,path] = uigetfile('*.mat','选择 WD＋VA 全局脊线实验结果');
    if isequal(name,0), fig=[]; return; end
    result = fullfile(path,name);
end
if ischar(result) || isstring(result)
    saved = load(result,'result');
    assert(isfield(saved,'result') && isfield(saved.result,'test'), ...
        '选择的文件不是 WD＋VA 全局脊线实验结果。');
    result = saved.result;
    loadedFromFile = true;
end
r = result.test;
if loadedFromFile || result.cfg.showFigures, visibility='on'; else, visibility='off'; end
fig = figure('Name','WD-RANSAC 与 WD＋VA 全局脊线基线', ...
    'Color','w','Visible',visibility,'Position',[80 80 1250 520]);
layout = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(layout);
plot(ax,r.snrDb,r.mseDb(:,1),'b-o','LineWidth',1.6); hold(ax,'on');
plot(ax,r.snrDb,r.mseDb(:,2),'r-s','LineWidth',1.6);
grid(ax,'on'); xlabel(ax,'SNR / dB（原 WD-RANSAC 口径）');
ylabel(ax,'最终轨迹 MSE / dB [(rad/时间单位)^2]');
title(ax,'最终直线轨迹误差'); legend(ax,r.methodNames,'Location','best');
ax = nexttile(layout);
plot(ax,r.snrDb,1000*r.meanTotalSeconds(:,1),'b-o','LineWidth',1.6); hold(ax,'on');
plot(ax,r.snrDb,1000*r.meanTotalSeconds(:,2),'r-s','LineWidth',1.6);
grid(ax,'on'); xlabel(ax,'SNR / dB（原 WD-RANSAC 口径）');
ylabel(ax,'单次平均总运行时间 / 毫秒');
title(ax,'包含 WD 和全部估计步骤'); legend(ax,r.methodNames,'Location','best');
title(layout,sprintf('F^{RR}、\\eta=0；固定 \\lambda=%.5g；M=%d；每个 SNR %d 次', ...
    r.lambda,result.cfg.mList,result.cfg.trials));
end
