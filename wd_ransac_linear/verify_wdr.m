function verify_wdr
% 数值检查入口：直接点击“运行”，全部通过后显示 PASS；不会写输出文件。
%VERIFY_WDR Numerical checks, runnable with the Editor's Run button.
% Checks mathematical conventions rather than matching one random plot.
scriptDir = fileparts(mfilename('fullpath')); % 本文件所在的 wd_ransac_linear 文件夹
moduleDir = fullfile(fileparts(scriptDir), 'modules'); % 旁边的 modules 文件夹
addpath(moduleDir);
cfg = wdr_config('quick');
rng(27,'twister');

% Positive/negative on-grid tones check frequency sign, factor 2, FFT shift.
for bin = [-100 -17 0 23 100]
    signal = wdr_generate_lfm(cfg,Inf,[0 bin*pi]);
    [wd,omega] = wdr_distribution(signal.x,cfg);
    ridge = wdr_ridge(wd,omega);
    assert(isequal(size(wd),[256 256]));
    assert(max(abs(ridge-bin*pi))<1e-10);
end

% Cross-check FFT WD against the explicit, symmetrized Eq.(1) sum.
signal = wdr_generate_lfm(cfg,-3,[220 -70]);
[wd,omega] = wdr_distribution(signal.x,cfg);
lag = (-128:128)';
weights = ones(size(lag)); weights([1 end]) = 0.5;
for col = [1 87 256]
    center = cfg.centerIndices(col);
    products = signal.x(center+lag).*conj(signal.x(center-lag));
    direct = exp(-1j*2*omega*(lag.'*cfg.dt))*(weights.*products);
    assert(max(abs(wd(:,col)-real(direct)))<1e-8);
    assert(max(abs(imag(direct)))<1e-8);
end

% A noiseless off-grid chirp has at most half a bin of ridge quantization.
clean = wdr_generate_lfm(cfg,Inf,[220 -70]);
[wd,omega,tWd] = wdr_distribution(clean.x,cfg);
ridge = wdr_ridge(wd,omega);
assert(max(abs(ridge-clean.omegaTrue))<=pi/2+1e-9);

% Noise convention: each real component has variance sigma^2.
nr = []; ni = [];
for k = 1:200
    noisy = wdr_generate_lfm(cfg,-5,[220 -70]);
    nr = [nr; real(noisy.noise)]; %#ok<AGROW>
    ni = [ni; imag(noisy.noise)]; %#ok<AGROW>
end
expectedVariance = 10^(5/10);
assert(abs(mean(nr.^2)/expectedVariance-1)<0.02);
assert(abs(mean(ni.^2)/expectedVariance-1)<0.02);

% Perfect ridge -> exact a,b and fully coherent amplitude N.
estimate = wdr_estimate(clean.x,cfg.t,tWd,clean.omegaTrue,[10 100],cfg);
assert(max(abs(estimate.a-220))<1e-9);
assert(max(abs(estimate.b+70))<1e-9);
assert(max(abs(estimate.score-cfg.nSamples))<1e-9);

% Deterministic pair/score reference uses all 512 samples, including ends.
pairs = [1 40 100; 256 210 180];
estimate = wdr_estimate(signal.x,cfg.t,tWd,ridge,[1 3],cfg,pairs);
manualScores = zeros(1,3);
for k = 1:3
    a = (ridge(pairs(2,k))-ridge(pairs(1,k))) / ...
        (tWd(pairs(2,k))-tWd(pairs(1,k)));
    b = ridge(pairs(1,k))-a*tWd(pairs(1,k));
    manualScores(k) = abs(sum(signal.x.*exp(-1j*(a*cfg.t.^2/2+b*cfg.t))));
end
assert(max(abs(manualScores-estimate.candidateScores))<1e-9);
assert(estimate.score(2)>=estimate.score(1));

% Sampling stays in bounds, excludes identical times, and is repeatable.
rng(901,'twister'); p1 = wdr_sample_pairs(256,3000,cfg);
rng(901,'twister'); p2 = wdr_sample_pairs(256,3000,cfg);
assert(isequal(p1,p2) && all(p1(1,:)<p1(2,:)));
assert(all(p1(:)>=1 & p1(:)<=256));

% End-to-end MC checks aggregation, repeatability, and parameter band.
cfg.trials = 3; cfg.snrDb = [-5 5]; cfg.mList = [10 100];
cfg.showProgress = false;
r1 = wdr_monte_carlo(cfg); r2 = wdr_monte_carlo(cfg);
assert(isequal(r1.ransacError,r2.ransacError));
assert(all(isfinite(r1.ransacMseDb(:))));
assert(all(abs(r1.trueB(:)-r1.trueA(:)/2)<=cfg.omegaLimit));
assert(all(abs(r1.trueB(:)+r1.trueA(:)/2)<=cfg.omegaLimit));
assert(abs(r1.ransacMse(1,1)-mean(r1.ransacError(:,1,1)))<1e-9);
fprintf('PASS: WD axis/sum, chirp ridge, noise, line fit, all-sample score, sampling and Monte Carlo.\n');
end
