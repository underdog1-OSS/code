function [aFine, bFine, detail] = wdr_oshea_refine(x, time, aCoarse, bCoarse, windowLength)
% 对线性调频候选做一次 O'Shea 精化，支持一批候选同时处理。
% 接入依据：QML-RANSAC 第 III.A 节 Step VI--X、式 (19)--(21)。
% 相位模型为 a*t^2/2+b*t；a 是角频率变化率，b 是 t=0 时的角频率。
% windowLength 是本实现可调的奇数平均长度，原文没有明确给出该数值。
% 本模块只接收粗估计和含噪信号，不读取真值。
x = x(:); time = time(:);
aCoarse = aCoarse(:).'; bCoarse = bCoarse(:).';
validateattributes(x,{'numeric'},{'vector','finite','nonempty'});
validateattributes(time,{'numeric'},{'vector','real','finite','numel',numel(x)});
validateattributes(aCoarse,{'numeric'},{'vector','real','finite','nonempty'});
validateattributes(bCoarse,{'numeric'},{'vector','real','finite','numel',numel(aCoarse)});
validateattributes(windowLength,{'numeric'},{'scalar','integer','positive','<=',numel(x)-2});
assert(mod(windowLength,2)==1,'移动平均长度必须是奇数。');
assert(all(diff(time)>0),'信号的时间轴必须严格递增。');
dt = mean(diff(time));
assert(max(abs(diff(time)-dt))<1e-10*max(1,abs(dt)), ...
    '居中移动平均要求信号等间隔采样。');

% 1. 用每条候选相位补偿原信号，使正确候选对应的残余频率接近零。
coarsePhase = (time.^2/2)*aCoarse+time*bCoarse;
dechirped = x.*exp(-1j*coarsePhase);

% 2. 沿时间方向做居中移动平均；仅使用完整窗口，不补零、不缩短端点窗。
% 每列是一条候选，卷积核只有一列，因此不会混合不同候选。
filtered = conv2(dechirped,ones(windowLength,1)/windowLength,'valid');
halfWindow = (windowLength-1)/2;
fitIndices = (1+halfWindow:numel(time)-halfWindow).';
fitTime = time(fitIndices);

% 3. 每列单独提取并展开残余相位。
residualPhase = unwrap(angle(filtered),[],1);

% 4. 拟合残余相位 d0+db*t+da*t^2/2，显式保留常数相位项。
% 用反斜杠解最小二乘，避免显式求逆；第三列的 1/2 与原相位模型一致。
basis = [ones(numel(fitTime),1),fitTime,fitTime.^2/2];
correction = basis\residualPhase;
bFine = bCoarse+correction(2,:);
aFine = aCoarse+correction(3,:);

% 粗估计不准时，相位展开和平均滤波都可能失效；是否采用由外层按评分决定。
if nargout>2
    detail.phase0 = correction(1,:);
    detail.deltaB = correction(2,:);
    detail.deltaA = correction(3,:);
    detail.fitIndices = fitIndices;
    detail.windowLength = windowLength;
end
end
