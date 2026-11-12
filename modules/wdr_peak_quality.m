function quality = wdr_peak_quality(wd, time, omega, peakIndex, options)
%WDR_PEAK_QUALITY 对每个 WD 峰点计算待验证的质量分数并排序。
% 输入 wd 为 [频率格 x 时刻]；peakIndex(i) 是第 i 列峰值的行号。
% z: 峰值相对背景的突出程度；c: 留一局部线性一致性；q = z*c。
% q 只安排抽样顺序，不是内点概率，也不参与相干评分和停止判断。
validateattributes(wd, {'numeric'}, {'real','finite','nonempty','2d'});
[nFreq,nTime] = size(wd);
validateattributes(time, {'numeric'}, {'real','finite','vector','numel',nTime});
validateattributes(omega, {'numeric'}, {'real','finite','vector','numel',nFreq});
validateattributes(peakIndex, {'numeric'}, ...
    {'real','finite','integer','vector','numel',nTime,'>=',1,'<=',nFreq});
time = time(:); omega = omega(:); peakIndex = peakIndex(:);
assert(all(diff(time)>0), 'wdr:QualityTime', ...
    'Time values must be strictly increasing and distinct.');
assert(all(diff(omega)>0), 'wdr:QualityFrequency', ...
    'Frequency values must be strictly increasing.');
validateattributes(options.guardBins, {'numeric'}, {'scalar','integer','nonnegative','finite'});
validateattributes(options.backgroundRadiusBins, {'numeric'}, {'scalar','integer','positive','finite'});
validateattributes(options.minBackgroundBins, {'numeric'}, {'scalar','integer','positive','finite'});
validateattributes(options.localNeighborsPerSide, {'numeric'}, {'scalar','integer','positive','finite'});
validateattributes(options.minLocalNeighbors, {'numeric'}, {'scalar','integer','>=',2,'finite'});
validateattributes(options.deltaFrequency, {'numeric'}, {'scalar','real','positive','finite'});
validateattributes(options.epsilonD, {'numeric'}, {'scalar','real','positive','finite'});
assert(options.backgroundRadiusBins>options.guardBins, 'wdr:QualityBackground', ...
    'Background radius must exceed the guard width.');
validateattributes(options.frequencyIndices, {'numeric'}, ...
    {'vector','integer','finite','nonempty','>=',1,'<=',nFreq});
search = options.frequencyIndices(:);
assert(all(diff(search)>0) && all(ismember(peakIndex,search)), ...
    'wdr:QualitySearch', 'Search indices must increase and contain every peak.');

quality.peakIndex = peakIndex;
quality.peakFrequency = omega(peakIndex);
quality.peakValue = zeros(nTime,1);
quality.backgroundMedian = zeros(nTime,1);
quality.backgroundScale = zeros(nTime,1);
quality.prominence = zeros(nTime,1);
quality.localResidual = nan(nTime,1);
quality.consistency = ones(nTime,1);
quality.localUsed = false(nTime,1);
quality.backgroundCount = zeros(nTime,1);
quality.backgroundFallback = false(nTime,1);

for i = 1:nTime
    k = peakIndex(i);
    distance = abs(search-k);
    background = search(distance>options.guardBins & distance<=options.backgroundRadiusBins);
    if numel(background)<options.minBackgroundBins
        background = search(distance>options.guardBins);
        quality.backgroundFallback(i) = true;
    end
    assert(numel(background)>=options.minBackgroundBins, 'wdr:QualityBackground', ...
        'Not enough background bins for time index %d.', i);
    values = wd(background,i);
    center = median(values);
    scale = 1.4826*median(abs(values-center));
    quality.peakValue(i) = wd(k,i);
    quality.backgroundMedian(i) = center;
    quality.backgroundScale(i) = scale;
    quality.backgroundCount(i) = numel(background);
    quality.prominence(i) = max(0,(wd(k,i)-center)/max(scale,options.epsilonD));

    % 留一：只用其他时刻的峰点拟合，避免被评价点参与自己的拟合。
    side = options.localNeighborsPerSide;
    left = max(1,i-side):i-1;
    right = i+1:min(nTime,i+side);
    neighbors = [left right];
    target = min(2*side,nTime-1);
    if numel(neighbors)<target
        remaining = setdiff([1:i-1 i+1:nTime],neighbors);
        [~,order] = sortrows([abs(time(remaining)-time(i)),remaining(:)],[1 2]);
        neighbors = [neighbors remaining(order(1:target-numel(neighbors)))]; %#ok<AGROW>
    end
    if numel(neighbors)>=options.minLocalNeighbors
        u = time(neighbors)-time(i);
        f = quality.peakFrequency(neighbors);
        slopes = zeros(numel(neighbors)*(numel(neighbors)-1)/2,1);
        p = 0;
        for j = 1:numel(neighbors)-1
            for k2 = j+1:numel(neighbors)
                p = p+1;
                slopes(p) = (f(j)-f(k2))/(u(j)-u(k2));
            end
        end
        localSlope = median(slopes);
        prediction = median(f-localSlope*u);
        residual = abs(quality.peakFrequency(i)-prediction);
        quality.localResidual(i) = residual;
        quality.consistency(i) = exp(-0.5*(residual/options.deltaFrequency)^2);
        quality.localUsed(i) = true;
    end
end
quality.q = quality.prominence.*quality.consistency;
assert(all(isfinite(quality.q)), 'wdr:QualityNonfinite', ...
    'Quality scores must be finite; check the WD scale and epsilonD.');
% 分数相同时按原时间编号排列，保证排序可复现。
[~,quality.order] = sortrows([-quality.q,(1:nTime)'],[1 2]);
quality.options = options;
end
