function [pairs, schedule] = wdr_prosac_pairs(rankOrder, count, growthBudget)
%WDR_PROSAC_PAIRS 按质量排序渐进抽取两个不同的时刻编号。
% rankOrder 从高分到低分存放原始时间编号。count 是固定抽样次数。
% growthBudget 是 PROSAC 的 T_N，控制前缀扩展日程，不是停止规则。
% 每轮抽样均在当前前缀中均匀进行；通过原时间编号返回给 wdr_estimate。
% 日程采用论文式 (3)-(6)：T_n=T_N*C(n,2)/C(N,2)，T'_2=1，
% T'_n=T'_(n-1)+ceil(T_n-T_(n-1))。日程内按式 (6) 包含新加入点；
% 超出完整前缀日程后，在所有点中均匀抽样。论文 Algorithm 1 对抽样
% 分支的叙述与式 (6) 有区别，此处明确采用式 (5)-(6) 的日程。
% 不使用原文基于内点概率的停止规则。
validateattributes(rankOrder, {'numeric'}, {'vector','integer','finite','nonempty'});
validateattributes(count, {'numeric'}, {'scalar','integer','positive','finite'});
validateattributes(growthBudget, {'numeric'}, {'scalar','integer','positive','finite'});
rankOrder = rankOrder(:);
nTime = numel(rankOrder);
assert(nTime>=2 && isequal(sort(rankOrder),(1:nTime)'), 'wdr:ProsacRanking', ...
    'Ranking must contain every time index exactly once and at least two points.');
threshold = zeros(nTime,1);
threshold(2) = 1;
for n = 3:nTime
    increment = growthBudget*2*(n-1)/(nTime*(nTime-1));
    threshold(n) = threshold(n-1)+ceil(increment);
end
pairs = zeros(2,count);
schedule.prefixSize = zeros(1,count);
schedule.forcedNewest = false(1,count);
schedule.threshold = threshold;
schedule.growthBudget = growthBudget;
schedule.rankOrder = rankOrder;
n = 2;
for t = 1:count
    while n<nTime && t>threshold(n)
        n = n+1;
    end
    schedule.prefixSize(t) = n;
    if t<=threshold(n)
        % 论文式 (6)：一个点来自前 n-1 个，另一个是前缀第 n 点。
        selectedRanks = [randperm(n-1,1),n];
        schedule.forcedNewest(t) = true;
    else
        % 完整前缀的扩展日程结束后，退化为完整集合内的均匀抽样。
        selectedRanks = randperm(n,2);
    end
    pairs(:,t) = sort(rankOrder(selectedRanks));
end
end
