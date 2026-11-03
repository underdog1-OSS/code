function [wd, omega, time] = wdr_distribution(x, cfg) % 输入信号x，参数配置结构体cfg，计算出的wd矩阵，wd每一行对应的频率omega，每一列对应的时间time 
% WD 模块：计算式 (1)，返回 256 个频率点 × 256 个时刻的时频矩阵。
% 注意指数中有 2*omega，因此频率轴间隔为 pi/(L*dt)，不要多乘一个 2。
%WDR_DISTRIBUTION Eq.(1), integer-lag discrete Wigner distribution.
% wd: [frequency bin x time index], NOT abs(WD) or abs(WD).^2.
% R(k,n)=x(n+k)*conj(x(n-k)), exponential is exp(-1j*2*omega*k*dt).
% Hence omega = pi*q/(L*dt), not 2*pi*q/(L*dt). //test
% For L=256, dt=1/256: omega=(-128:127)*pi, bin spacing=pi.
% The even-length lag endpoint is averaged with its conjugate counterpart:
% this equals half weights on k=+/-128 and full weights on k=-127:127.
% Only the fully supported middle interval is evaluated; no edge padding.
x = x(:); % 把输入排成一列
assert(numel(x)==cfg.nSamples, 'wdr:SignalSize', 'Unexpected signal length.'); % 检查x是否有512个元素
L = cfg.wdLength;   % 取出配置的WD长度（256） 
assert(mod(L,2)==0, 'wdr:WdLength', 'WD length must be even.');  % 检查L是否为偶数
lag = (-L/2:L/2-1)'; % 相对中心采样点的偏移编号
center = cfg.centerIndices(:).'; % 取出中间256个点的位置编号
plus = bsxfun(@plus, lag, center); % lag和center交叉相加，得到256*256矩阵，每一列表示同一个中心点，每一行表示同一个偏移值
minus = bsxfun(@minus, center, lag); % center和lag交叉相减，得到256*256矩阵
assert(min(plus(:))>=1 && max(plus(:))<=numel(x) && ... % 检查plus和minus中的所有编号在1-512之间
    min(minus(:))>=1 && max(minus(:))<=numel(x), ...
    'wdr:WindowBounds', 'Lag window exceeds available signal samples.');
correlation = x(plus).*conj(x(minus)); % 取中心点右侧对应位置的信号值，乘中心点左侧对应位置的信号值的共轭相当于论文中的x（t+kΔt）* x*（t-kΔt）
correlation(1,:) = real(correlation(1,:));
wd = real(fftshift(fft(ifftshift(correlation,1),L,1),1));% 计算wd
omega = (-L/2:L/2-1)'*pi/(L*cfg.dt); % 生成频率坐标
time = cfg.t(cfg.centerIndices);     % 生成时间坐标
end

