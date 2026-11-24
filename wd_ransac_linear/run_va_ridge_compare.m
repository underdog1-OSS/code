%% WD＋VA 全局脊线＋普通最小二乘，与原 WD-RANSAC 的独立配对对照。
% 先运行数学与无噪声验证，再独立训练固定 lambda，最后独立测试。
% 默认只显示一张包含两个子图的大图，数据留在工作区 result，不自动保存文件。
scriptDir = fileparts(mfilename('fullpath'));
codeDir = fileparts(scriptDir);
addpath(fullfile(codeDir,'modules'));
cfg = wdr_va_config;             % 试验次数等参数集中在这个配置模块。
verification = wdr_va_verify(cfg);
result = wdr_va_experiment(cfg);
result.verification = verification;

% 将实际运行所用的代码文本放进 MAT，避免另外生成一套源代码文件夹。
sourceNames = {'wdr_va_config','wdr_va_ridge','wdr_va_verify','wdr_va_experiment', ...
    'wdr_va_plot','wdr_config','wdr_distribution','wdr_generate_lfm', ...
    'wdr_estimate','wdr_ridge','wdr_sample_pairs'};
for k = 1:numel(sourceNames)
    result.source(k).name = [sourceNames{k} '.m'];
    result.source(k).text = fileread(fullfile(codeDir,'modules',[sourceNames{k} '.m']));
end
result.entryScript = fileread([mfilename('fullpath') '.m']);
% 如主动开启数值保存，则手动选择一个 MAT 文件，不创建文件夹。
if cfg.saveResults
    [dataFile,dataPath] = uiputfile('*.mat','选择数值结果保存位置', ...
        'WD_VA_GlobalRidge_results.mat');
    if isequal(dataFile,0)
        fprintf('已跳过数值保存，结果仍保留在工作区 result。\n');
    else
        save(fullfile(dataPath,dataFile),'result','-v7.3');
        fprintf('数值结果已保存：%s\n',fullfile(dataPath,dataFile));
    end
end
fig = wdr_va_plot(result);
disp(result.test.summary);
fprintf('\n实验完成：合并图已显示，数据保留在工作区 result，图片请手动保存。\n');
fprintf('若已手动保存 MAT，可打开 wdr_va_plot.m 点击运行并选择该文件重绘。\n');
