function verification = wdr_va_verify(cfg)
% 核对全局递推、全起点初始化、回溯、频率单位及线性模型的无噪声结果。
if nargin<1, cfg = wdr_va_config; end
state = rng; cleaner = onCleanup(@() rng(state)); %#ok<NASGU>
rng(cfg.seed+7000000,'twister');
messages = {};
toy = [4 1 3 1;2 5 1 2;1 2 4 6];
for lambda = [0 0.2 1]
    va = wdr_va_ridge(toy,(0:2)',lambda,cfg.vaFloorRelative);
    % 穷举所有短路径，不使用 VA 的前驱信息作为参考答案。
    for n = 1:4
        bestByEnd = inf(3,1);
        for code = 0:3^n-1
            digits = code; path = zeros(n,1);
            for k = 1:n
                path(k) = mod(digits,3)+1; digits = floor(digits/3);
            end
            cost = sum(va.nodeCost(sub2ind(size(toy),path,(1:n).')))+ ...
                lambda*sum(diff(path).^2);
            bestByEnd(path(end)) = min(bestByEnd(path(end)),cost);
        end
        assert(max(abs(bestByEnd-va.cumulativeCost(:,n)))<1e-10);
    end
    tracedCost = va.dataCost+va.transitionCost;
    assert(abs(tracedCost-va.totalCost)<1e-9);
end
messages{end+1} = '通过：所有时间层、所有终点的累计代价与穷举一致，回溯代价正确。';
% 全局最优从首列的次高峰起步，验证没有强制跟随首列最大峰。
va = wdr_va_ridge([10 1e-4 1e-4;1 100 100],[0;1],3);
assert(all(va.index==2));
messages{end+1} = '通过：保留全部起点，能选出首列不是最大峰的全局最优路径。';
v0 = wdr_va_ridge(toy,(0:2)',0);
[~,maximum] = max(toy,[],1);
assert(isequal(v0.index,maximum(:)));
v1 = wdr_va_ridge(toy,(0:2)',0.2);
v2 = wdr_va_ridge(1e6*toy,(0:2)',0.2);
assert(isequal(v1.index,v2.index));
zero = wdr_va_ridge(zeros(4,5),(0:3)',0.1);
negative = wdr_va_ridge([-1 0;2 -3;0 4],(0:2)',0.1);
assert(all(isfinite(zero.nodeCost(:))) && all(zero.index==1));
assert(all(isfinite(negative.nodeCost(:))));
messages{end+1} = '通过：lambda=0 退化为最大峰，整体幅度缩放不改变路径，负值及全零不产生 log(0)。';

signal = wdr_generate_lfm(cfg,Inf,cfg.paperParameters);
[wd,omega,tWd] = wdr_distribution(signal.x,cfg);
assert(isequal(size(wd),[256 256]));
assert(max(abs(diff(omega)-pi))<1e-10);
assert(min(cfg.centerIndices)-128>=1 && max(cfg.centerIndices)+128<=numel(signal.x));
va = wdr_va_ridge(wd,omega,0.032,cfg.vaFloorRelative);
verification.noiselessRawMaxError = max(abs(va.omega-signal.omegaTrue));
assert(verification.noiselessRawMaxError<=pi/2+1e-8);
coefficient = [tWd ones(size(tWd))]\va.omega;
fit = [tWd ones(size(tWd))]*coefficient;
verification.noiselessFitMse = mean((fit-signal.omegaTrue).^2);
assert(verification.noiselessFitMse<=mean((va.omega-signal.omegaTrue).^2)+1e-8);
% 改变实际频率单位不应改变格点路径；lambda 的物理尺度要相应换算。
hz = wdr_va_ridge(wd,omega/(2*pi*cfg.secondsPerTimeUnit),0.032,cfg.vaFloorRelative);
assert(isequal(va.index,hz.index));
assert(abs(hz.lambdaOmega-va.lambdaOmega*(2*pi*cfg.secondsPerTimeUnit)^2)<1e-9);
messages{end+1} = '通过：论文线性信号的 WD 轴与完整边界正确，无噪声原始脊线误差不超过半格。';
high = wdr_generate_lfm(cfg,30,cfg.paperParameters);
[highWd,highOmega,highTime] = wdr_distribution(high.x,cfg);
highPath = wdr_va_ridge(highWd,highOmega,0.032,cfg.vaFloorRelative);
highDesign = [highTime ones(size(highTime))];
highFit = highDesign*(highDesign\highPath.omega);
verification.highSnrFitMse = mean((highFit-high.omegaTrue).^2);
assert(verification.highSnrFitMse<pi^2,'高 SNR 下未恢复到一个频率格点的精度。');
messages{end+1} = '通过：30 dB 线性信号的拟合轨迹达到一个频率格点以内的均方精度。';
verification.messages = messages;
verification.passed = true;
verification.paperLinearParameters = cfg.paperParameters;
verification.timeRange = tWd([1 end]);
verification.frequencyRange = omega([1 end]);
verification.deltaOmega = pi;
verification.secondsPerTimeUnit = cfg.secondsPerTimeUnit;
for k = 1:numel(messages), fprintf('%s\n',messages{k}); end
end
