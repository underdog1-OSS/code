function cfg = wdr_config(preset) % 参数配置 
% 配置模块：采样间隔、时间范围、WD 长度、SNR、迭代次数等集中在这里。
% 原文设定：512 点信号、中心 256 点 WD、256 点窗口、每个 SNR 1000 次。
% 请注意：论文 SNR 的分母是“实部或虚部各自的方差”，不是复噪声总功率。
% 原文未给出的细节，本实现采用以下明确约定：
%   1) SNR 间隔 1 dB；随机生成左右端点频率，均匀覆盖允许的参数区域。
%   2) 高斯抽点均值为中心索引，标准差为时刻数/4，越界或重合时重新抽取。
%   3) WD 使用矩形延迟窗，偶数长度端点采用正负延迟各一半的对称处理。
%   4) 不同 M 共用同一随机序列的前 M 对，便于公平比较。
% 因而这是按论文描述实现的复现代码，不能保证与原作者曲线逐点重合。
% MSE 对中心 256 点角频率误差和所有试验求平均；单位不是 Hz^2。
% WD 暂用周期误差作对照实验；WDRANSAC 使用直接误差，论文定义仍需核实。
% 原文斜率分母 n2-n2 是排印问题，本实现使用两点斜率所需的 n2-n1。
%WDR_CONFIG Settings for Djurovic (2016), Example 1, MATLAB R2020b.
% Source: A WD-RANSAC Instantaneous Frequency Estimator,
% IEEE Signal Processing Letters 23(5):757-761, Sec. II-III and Example 1.
% DOI: 10.1109/LSP.2016.2551732.
%
% PAPER-SPECIFIED:
%   x(t) = exp(1j*(a*t^2/2+b*t)); t in [-1,1), dt=1/256, N=512.
%   WD at t in [-1/2,1/2), 256 lag samples, a 256-by-256 TF matrix.55
%   |b-a/2| and |b+a/2| <= pi/(2*dt).
%   Noise real/imaginary parts INDEPENDENT, each with variance sigma^2.
%   Paper SNR = 10*log10(A^2/sigma^2), not A^2/(2*sigma^2).
%   1000 trials/SNR. Fig. 1: M=10,100,300,900,3000.
%   Two different Gaussian-distributed instants per candidate.
%   Score = abs(sum(x.*exp(-1j*(aHat*t.^2/2+bHat*t)))) over ALL 512 samples.
%
% EXPLICIT IMPLEMENTATION CHOICES (not supplied by the paper):
%   (1) SNR spacing = 1 dB (paper specifies only range [-15,5]).
%   (2) Draw the IF at t=-0.5 and +0.5 independently uniformly in the
%       admissible band. This is uniform over the admissible (a,b) diamond.
%   (3) Gaussian time indices: mean=(K+1)/2, std=K/4; round to integer,
%       reject out-of-range/repeated indices. No minimum pair separation.
%   (4) Rectangular lag window; at the even FFT's unique Nyquist lag,
%       average the +/-L/2 products to preserve a real-valued WD.
%   (5) All M budgets use prefixes of the same candidate sequence per trial.
%       This couples comparisons but keeps each budget's estimator unchanged.
%   (6) MSE = mean squared IF error on the middle 256 samples, averaged
%       over trials BEFORE 10*log10. Errors use rad/time-unit, not Hz.
%
% The printed two-point slope denominator on p.758 has n2-n2; the code
% uses n2-n1 as required by the line equation. There is no LS refit,
% inlier threshold, candidate-feasibility filter, or local optimization.
if nargin < 1, preset = 'paper'; end % 选择paper模式，按照论文的模式进行仿真，还有quick模式是减少蒙特卡洛次数来进行测试
cfg.preset = validatestring(preset, {'paper','quick'});
cfg.dt = 1/256;              % 采样间隔
cfg.nSamples = 512;          % 信号的512个点 
cfg.wdLength = 256;          % 计算WD使用256个点
cfg.t = (-256:255)' * cfg.dt;% 每个采样点对应的时间
cfg.centerIndices = (129:384)';% 中间256个点的位置编号
cfg.tWd = cfg.t(cfg.centerIndices);% 中间256个采样点对应的时间
cfg.omegaLimit = pi/(2*cfg.dt);% k每增加1，相位变化2*omega*dt，要保证在-pi到pi之间，算出omega的边界
cfg.amplitude = 1;            % 复指数信号幅度为1
cfg.phase0 = 0;               % 信号额外的初始相位为0
cfg.snrDb = -15:1:5;          % 论文定义下的 SNR 范围
cfg.trials = 1000;            % 每个 SNR 的蒙特卡洛次数
cfg.mList = [10 100 300 900 3000]; % 图 1 比较的随机点对数量（轮数）
cfg.seed = 20261023;          % 随机数种子
cfg.gaussianMeanIndex = (numel(cfg.centerIndices)+1)/2; % 中间共有256个点，选择选点倾向的中心，numel计算向量的元素数
cfg.gaussianStdIndex = numel(cfg.centerIndices)/4;      % 选点散开的程度
cfg.scoreBlockSize = 256;     % 计算评分时每批最多处理256条候选直线
cfg.showProgress = true;      % 一个开关，表示运行多次实验时显示进度，改成False不显示
if strcmp(cfg.preset, 'quick')% 选择quick模式，修改蒙特卡洛次数和信噪比范围
    cfg.trials = 50;
    cfg.snrDb = -15:2:5;
end
end

