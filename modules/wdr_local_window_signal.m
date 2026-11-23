function signal = wdr_local_window_signal(cfg, snrDb, physicalParameters)
% 复用原信号生成器；把物理起点的 f0、斜率转换到原来的居中时间轴。
% 输入参数为 [f0(Hz), a(Hz/s)]，省略时按本实验范围均匀生成。
if nargin < 3 || isempty(physicalParameters)
    f0 = cfg.startFrequencyRangeHz(1)+diff(cfg.startFrequencyRangeHz)*rand;
    rate = cfg.chirpRateRangeHzPerSec(1)+diff(cfg.chirpRateRangeHzPerSec)*rand;
else
    f0 = physicalParameters(1);
    rate = physicalParameters(2);
end
origin = cfg.t(1);
angularParameters = 2*pi*[rate, f0-rate*origin];
signalCfg = cfg;
% 令物理时间零点的相位等于 cfg.phase0，补偿换时间原点产生的常相位。
signalCfg.phase0 = cfg.phase0+pi*rate*origin^2-2*pi*f0*origin;
signal = wdr_generate_lfm(signalCfg,snrDb,angularParameters);
signal.f0Hz = f0;
signal.rateHzPerSec = rate;
signal.physicalTime = cfg.t-origin;
signal.trueHz = signal.omegaTrue/(2*pi);
assert(all(f0+rate*signal.physicalTime>0), '真实信号必须始终为正频率。');
assert(max(f0+rate*signal.physicalTime)<cfg.omegaLimit/(2*pi)-2, ...
    '真实信号距 WD 频带边缘必须至少留出 2 Hz。');
end
