function [estimate,va,schedule,timing] = wdr_va_prosac_estimate(x,time,tWd,wd,omega,mList,cfg,lambda,epsilonP,growthBudget)
% 独立入口：VA 选点 -> PROSAC 时间点对 -> 原相干 J。
timer=tic; va=wdr_va_select(wd,omega,tWd,lambda,epsilonP); timing.vaSeconds=toc(timer);
timer=tic; [pairs,schedule]=wdr_prosac_pairs(va.rankOrder,max(mList),growthBudget); timing.samplingSeconds=toc(timer);
timer=tic; estimate=wdr_estimate(x,time,tWd,va.selectedOmega,mList,cfg,pairs); timing.scoringSeconds=toc(timer);
schedule.pairs=pairs;
end
