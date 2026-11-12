%% 比较不同候选数量 M 下，三种抽样方法的瞬时频率误差
% 在 MATLAB 编辑器中点击“运行”，不自动保存图像或数据文件。
% 建议图像名称：WD_RANSAC_PROSAC_M_comparison_50trials_TN3000.svg
% 如果已有 run_prosac_compare.m 得到的 result，则直接使用，不重新仿真。
% 如果没有可用的 result，则自动运行 quick 模式，每个 SNR 50 次。
% 五个 M 共用同一批候选及同一 PROSAC 增长日程 T_N，不是五次独立实验。
scriptDir = fileparts(mfilename('fullpath'));
moduleDir = fullfile(fileparts(scriptDir),'modules');
addpath(moduleDir);

requiredFields = {'cfg','growthBudget','snrDb','mList', ...
    'gaussianMseDb','uniformMseDb','prosacMseDb'};
requestedM = [10 100 300 900 3000];
useExisting = exist('result','var')==1 && isstruct(result) && ...
    isscalar(result) && all(isfield(result,requiredFields));
if useExisting
    useExisting = isstruct(result.cfg) && isfield(result.cfg,'trials') && ...
        all(ismember(requestedM,result.mList)) && ...
        isequal(size(result.gaussianMseDb),[numel(result.snrDb),numel(result.mList)]) && ...
        isequal(size(result.uniformMseDb),size(result.gaussianMseDb)) && ...
        isequal(size(result.prosacMseDb),size(result.gaussianMseDb));
end
if useExisting
    mCompareData = result;
    fprintf('Using the existing result in the workspace; no new simulation.\n');
else
    fprintf('No compatible result in workspace. Running the quick experiment (50 trials/SNR).\n');
    mCompareCfg = wdr_config('quick');
    mCompareData = wdr_compare_prosac(mCompareCfg,max(mCompareCfg.mList));
end

% 用统一纵轴范围，方便比较不同 M；均先对线性 MSE 求平均再转换 dB。
[~,mColumns] = ismember(requestedM,mCompareData.mList);
allDb = [mCompareData.gaussianMseDb(:,mColumns), ...
         mCompareData.uniformMseDb(:,mColumns), ...
         mCompareData.prosacMseDb(:,mColumns)];
finiteDb = allDb(isfinite(allDb));
assert(~isempty(finiteDb),'wdr:PlotNoFiniteData','No finite MSE values to plot.');
yLimits = [floor(min(finiteDb)/10)*10-5, ceil(max(finiteDb)/10)*10+5];

figure('Name','WD-RANSAC sampling comparison across M','Color','w', ...
    'Position',[100 100 1450 800]);
for k = 1:numel(requestedM)
    col = mColumns(k);
    subplot(2,3,k);
    plot(mCompareData.snrDb,mCompareData.gaussianMseDb(:,col),'b-o','LineWidth',1.2); hold on;
    plot(mCompareData.snrDb,mCompareData.uniformMseDb(:,col),'k-s','LineWidth',1.2);
    plot(mCompareData.snrDb,mCompareData.prosacMseDb(:,col),'r-^','LineWidth',1.2);
    grid on; xlabel('SNR (dB, paper convention)'); ylabel('Direct IF MSE (dB)');
    ylim(yLimits);
    title(sprintf('M = %d',requestedM(k)));
end

% 第六张小图：固定一个 SNR，横向比较不同 M。若没有 -7 dB，取最近值。
targetSnr = -7;
[~,snrIndex] = min(abs(mCompareData.snrDb-targetSnr));
subplot(2,3,6);
semilogx(requestedM,mCompareData.gaussianMseDb(snrIndex,mColumns),'b-o','LineWidth',1.2); hold on;
semilogx(requestedM,mCompareData.uniformMseDb(snrIndex,mColumns),'k-s','LineWidth',1.2);
semilogx(requestedM,mCompareData.prosacMseDb(snrIndex,mColumns),'r-^','LineWidth',1.2);
grid on; xlabel('Number of candidates M'); ylabel('Direct IF MSE (dB)');
xticks(requestedM); ylim(yLimits);
title(sprintf('SNR = %g dB',mCompareData.snrDb(snrIndex)));
legend('Original Gaussian','Uniform','Quality + PROSAC','Location','best');
sgtitle(sprintf('Fixed sampling budgets: %d trials/SNR, T_N = %d', ...
    mCompareData.cfg.trials,mCompareData.growthBudget));
fprintf('Figure displays M = 10, 100, 300, 900, 3000; each uses the first M candidates.\n');
fprintf('Suggested filename: WD_RANSAC_PROSAC_M_comparison_%dtrials_TN%d.svg\n', ...
    mCompareData.cfg.trials,mCompareData.growthBudget);
