function pairs = wdr_sample_pairs(nTime, count, cfg)                                % nTime：ridge 上可选的时间点数量 = 256
% 抽样模块：每轮抽两个不同的时刻索引，采用高斯抽样并排序，属于WDRANSAC中选点求候选直线  % count：要抽多少对编号，例如要生成 3000 条候选直线，就抽 3000 对
% WDR_SAMPLE_PAIRS Step 3: two distinct rounded Gaussian time indices.      
% Out-of-range draws are rejected, not clipped to the interval endpoints.
% Repetition of the same pair in different iterations is allowed.
validateattributes(count, {'numeric'}, {'scalar','integer','positive'});
assert(nTime>=2 && cfg.gaussianStdIndex>0, 'wdr:PairSettings', ...
    'At least two time instants and positive Gaussian std are required.');
assert(cfg.gaussianMeanIndex>=1 && cfg.gaussianMeanIndex<=nTime, ...
    'wdr:PairMean', 'Gaussian mean must be inside the time-index interval.');
pairs = zeros(2,count);                                             % 2*count的矩阵，表示选择的count对点
pending = 1:count;                                                  % 记录那些列还没选好
while ~isempty(pending)                                             % 只要有没选好的列就循环
    draw = round(cfg.gaussianMeanIndex + ...                        % 为每个待填的列生成两个随机数；乘标准差、加中心值后，编号会更倾向于落在高斯分布中心附近
        cfg.gaussianStdIndex*randn(2,numel(pending)));              % draw是本轮抽到的候选编号（2*n矩阵，n是成功抽出的列数）
    good = all(draw>=1 & draw<=nTime,1) & draw(1,:)~=draw(2,:);
    pairs(:,pending(good)) = sort(draw(:,good),1);
    pending = pending(~good);
end
end

