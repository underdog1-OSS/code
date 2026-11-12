function options = wdr_quality_config(wd, omega)
%WDR_QUALITY_CONFIG 质量排序的初版参数，可在主脚本中修改。
% WD 的行对应频率、列对应时间；omega 的单位是角频率。
% 这些数值是卡片中的实验初值，并非论文已经验证的最优值。
validateattributes(wd, {'numeric'}, {'real','finite','nonempty','2d'});
validateattributes(omega, {'numeric'}, {'real','finite','vector','numel',size(wd,1)});
omega = omega(:);
assert(numel(omega)>=2 && all(diff(omega)>0), 'wdr:QualityFrequency', ...
    'Frequency grid must contain at least two strictly increasing values.');
options.guardBins = 2;                 % 避开主峰及其左右相邻 2 格
options.backgroundRadiusBins = 16;    % 在主峰两侧最多 16 格内找背景
options.minBackgroundBins = 8;
options.localNeighborsPerSide = 3;
options.minLocalNeighbors = 4;
options.deltaFrequency = 2*median(diff(omega)); % 残差尺度：初始为两个频率格
% 根据当前 WD 的数值量级设置下限；全零 WD 时仍保证下限为正。
options.epsilonD = max(1e-12*max(abs(wd(:))),realmin('double'));
options.frequencyIndices = (1:numel(omega))'; % 与原算法一致：搜索全部频率格
end
