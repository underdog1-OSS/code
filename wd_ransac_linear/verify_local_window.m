function verification = verify_local_window
% 数值检查：无噪声坐标、边界、窗口递推、共用抽样和原相干评分。
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(scriptDir),'modules'));
cfg = wdr_local_window_config;
rngState = rng;
cleaner = onCleanup(@() rng(rngState)); %#ok<NASGU>
rng(20261123,'twister');
messages = {};
signal = wdr_local_window_signal(cfg,Inf,[20.25 8]);
[wd,omega,tWd] = wdr_distribution(signal.x,cfg);
assert(isequal(size(wd),[256 256]));
assert(max(abs(diff(omega)/(2*pi)-0.5))<1e-12);
assert(abs(tWd(1)-cfg.t(1)-0.5)<1e-12);
assert(abs(tWd(end)-cfg.t(1)-1.49609375)<1e-12);
assert(min(cfg.centerIndices)-cfg.wdLength/2>=1);
assert(max(cfg.centerIndices)+cfg.wdLength/2<=cfg.nSamples);
expected = cfg.amplitude*exp(1j*(2*pi*20.25*cfg.physicalTime+ ...
    pi*8*cfg.physicalTime.^2+cfg.phase0));
assert(max(abs(expected-signal.clean))<1e-10);
messages{end+1} = '通过：512 点信号、256×256 WD、0.5 Hz 间隔和完整延迟支持区间。';
verification.maxNoiselessPeakErrorHz = zeros(3,1);
for m = 1:3
    ridge = wdr_local_window_ridge(wd,omega,cfg,char('A'+m-1),signal.omegaTrue(1));
    verification.maxNoiselessPeakErrorHz(m) = max(abs(ridge/(2*pi)-signal.trueHz));
    assert(verification.maxNoiselessPeakErrorHz(m)<=cfg.frequencyStepHz/2+1e-9);
end
messages{end+1} = '通过：无噪声三组峰值偏差均不超过半个格点，即 0.25 Hz。';
% 独立求和核对首、中、末有效时间列的 WD，不仅检查画出的坐标。
lag = (-128:128)'; weights = ones(size(lag)); weights([1 end]) = 0.5;
for col = [1 128 256]
    center = cfg.centerIndices(col);
    products = signal.x(center+lag).*conj(signal.x(center-lag));
    direct = exp(-1j*2*omega*(lag.'*cfg.dt))*(weights.*products);
    assert(max(abs(real(direct)-wd(:,col)))<1e-8);
    assert(max(abs(imag(direct)))<1e-8);
end
messages{end+1} = '通过：有效区间首、中、末列与直接 WD 加权求和一致。';
% 手造三个截面：全频最大峰在 10 Hz，局部峰应由 5 递推到 5.5、6 Hz。
grid = (0:0.5:12)'*2*pi;
toy = zeros(numel(grid),3);
toy(11,1)=20; toy(21,2:3)=100; toy(12,2)=10; toy(13,3)=10;
[b,detail] = wdr_local_window_ridge(toy,grid,cfg,'B');
assert(max(abs(b/(2*pi)-[5;5.5;6]))<1e-10);
assert(abs(detail.lower(3)/(2*pi)-4.95)<1e-10);
assert(abs(detail.upper(3)/(2*pi)-6.05)<1e-10);
messages{end+1} = '通过：局部中心依次为 5、5.5 Hz，不会被后续 10 Hz 全频噪声峰重置。';
% 首列误选到 2 Hz 后，B 不会偷偷用真值纠错；C 仅改首列初始化。
toy = zeros(numel(grid),12); toy(5,:)=10; toy(11,:)=100; toy(5,1)=200;
b = wdr_local_window_ridge(toy,grid,cfg,'B');
c = wdr_local_window_ridge(toy,grid,cfg,'C',5*2*pi);
assert(all(abs(b/(2*pi)-2)<1e-10) && all(abs(c/(2*pi)-5)<1e-10));
messages{end+1} = '通过：错误初始化产生持续错误段；理想初始化组只在首列读取真值。';
% 零频窗口退化和负频初始化也按原规则处理，不删除试验。
toy(1,1)=300;
b = wdr_local_window_ridge(toy,grid,cfg,'B');
assert(all(b==0));
negativeGrid = (-12:0.5:12)'*2*pi;
toy = zeros(numel(negativeGrid),2); toy(15,1)=10; toy(14,2)=10;
b = wdr_local_window_ridge(toy,negativeGrid,cfg,'B');
assert(max(abs(b/(2*pi)-[-5;-5.5]))<1e-10);
messages{end+1} = '通过：零频单格锁定、负频窗口上下界排序，均无回退或重启。';
% A 的峰和旧函数逐点一致，仍然取有符号 WD 最大值。
a = wdr_local_window_ridge(wd,omega,cfg,'A');
reference = wdr_ridge(wd,omega);
assert(isequal(a,reference));
noisy = wdr_local_window_signal(cfg,-7,[20.25 8]);
[noisyWd,omega] = wdr_distribution(noisy.x,cfg);
a = wdr_local_window_ridge(noisyWd,omega,cfg,'A');
assert(isequal(a,wdr_ridge(noisyWd,omega)));
messages{end+1} = '通过：A 与现有全频最大峰模块逐点相同。';
pairs = [1 40 100;256 210 180];
estimate = wdr_estimate(noisy.x,cfg.t,tWd,a,3,cfg,pairs);
manualScores = zeros(1,3);
for k = 1:3
    i=pairs(1,k); j=pairs(2,k);
    slope=(a(j)-a(i))/(tWd(j)-tWd(i));
    intercept=a(i)-slope*tWd(i);
    manualScores(k)=abs(sum(noisy.x.*exp(-1j*(slope*cfg.t.^2/2+intercept*cfg.t))));
end
assert(max(abs(manualScores-estimate.candidateScores))<1e-9);
perfect = wdr_estimate(signal.x,cfg.t,tWd,signal.omegaTrue,3,cfg,pairs);
assert(abs(perfect.a/(2*pi)-8)<1e-9);
assert(abs((perfect.b+perfect.a*cfg.t(1))/(2*pi)-20.25)<1e-9);
assert(abs(perfect.score-cfg.nSamples*cfg.amplitude)<1e-8);
messages{end+1} = '通过：相干评分使用全部 512 点，物理 f0、斜率转换及理想直线评分正确。';
% 小型端到端复现，验证指标定义、种子、共用抽样和表格汇总。
small = cfg; small.trials=2; small.snrDb=-7; small.mList=20;
small.bootstrapCount=50; small.showProgress=false;
r1=wdr_local_window_compare(small); r2=wdr_local_window_compare(small);
assert(isequal(r1.peakHz,r2.peakHz) && isequal(r1.trackMseHz2,r2.trackMseHz2));
assert(isequal(r1.initialCorrect(:,:,1),r1.initialCorrect(:,:,2)));
assert(all(r1.initialCorrect(:,:,3)==1));
assert(height(r1.summary)==3 && height(r1.trials)==6);
assert(abs(r1.summary.TrackMSE_Hz2(1)-mean(r1.trackMseHz2(:,1,1)))<1e-10);
assert(all(r1.totalSeconds(:)>0));
messages{end+1} = '通过：小型端到端结果可复现，初始化、计时和汇总表一致。';
verification.messages = messages;
verification.passed = true;
verification.validTimeSeconds = [0.5 1.49609375];
verification.frequencyBandHz = omega([1 end])/(2*pi);
verification.toleranceHz = cfg.correctToleranceHz;
for k=1:numel(messages), fprintf('%s\n',messages{k}); end
end
