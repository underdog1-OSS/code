function signal = wdr_generate_lfm(cfg, snrDb, parameters)
% 信号模块：生成 x(t)=exp(j*(a*t^2/2+b*t))，并按论文 SNR 加复高斯噪声。
%WDR_GENERATE_LFM Generate an LFM signal with the paper's noise convention.
% parameters = [a b] is optional; omitted -> random admissible parameters.
% snrDb=Inf produces a noiseless signal for numerical checks.
if nargin < 3 || isempty(parameters)                                        % 没传parameters
    endpoints = (2*rand(2,1)-1)*cfg.omegaLimit;                             % 这一行先随机生成两个角频率：分别作为两个边界处的值。
    a = endpoints(2)-endpoints(1); % Endpoints at -0.5 and +0.5, span=1     % 计算这个时间-频率点对应的a，b
    b = mean(endpoints);
else
    validateattributes(parameters, {'numeric'}, {'real','finite','numel',2});
    a = parameters(1); b = parameters(2);
end
assert(max(abs([b-a/2 b+a/2])) <= cfg.omegaLimit + 1e-10, ...               % 检查中间 WD 时间范围两端的真实频率，是否都在允许范围内
    'wdr:AliasBand', 'True IF must satisfy |b +/- a/2| <= pi/(2*dt).');
signal.t = cfg.t;                                                           % 保存时间t和真实的a，b
signal.a = a; signal.b = b;
signal.clean = cfg.amplitude*exp(1j*(a*cfg.t.^2/2+b*cfg.t+cfg.phase0));     % 生成未加噪声的线性调频信号
signal.sigma = cfg.amplitude*10^(-snrDb/20);                                % 计算噪声标准差
signal.noise = signal.sigma*(randn(size(cfg.t))+1j*randn(size(cfg.t)));     % 生成复高斯噪声
signal.x = signal.clean+signal.noise;                                       % 生成带噪信号
signal.omegaTrue = a*cfg.tWd+b;                                             % 计算真实角频率
signal.snrDb = snrDb;                                                       % 保存信噪比
end

