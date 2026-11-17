function va = wdr_va_select(wd, omega, time, lambda, epsilonP)
%WDR_VA_SELECT 从每一列的 VA 最小累计代价选点并按平均代价排序。
% wd 为频率 x 时间；不回溯最终路径，也不更改原 WD 最大峰。
validateattributes(wd,{'numeric'},{'real','finite','nonempty','2d'});
[nFreq,nTime] = size(wd);
validateattributes(omega,{'numeric'},{'real','finite','vector','numel',nFreq});
validateattributes(time,{'numeric'},{'real','finite','vector','numel',nTime});
validateattributes(epsilonP,{'numeric'},{'real','finite','scalar','positive','<=',1});
omega=omega(:); time=time(:);
assert(nFreq>=2 && all(diff(omega)>0), 'wdr:VaFrequency', 'Frequency grid must increase.');
spacing=diff(omega);
assert(max(abs(spacing-spacing(1)))<=1e-9*abs(spacing(1)), ...
    'wdr:VaFrequency', 'Frequency grid must be uniformly spaced.');
assert(all(diff(time)>0), 'wdr:VaTime', 'Time must strictly increase.');
p=max(wd,0); total=sum(p(:));
assert(isfinite(total) && total>0, 'wdr:VaNonpositive', ...
    'WD must contain finite positive energy for global normalization.');
normalized=p/total;
va.nodeCost=-log(max(normalized,epsilonP));
[va.cumulativeCost,va.predecessor]=wdr_va_recurrence(va.nodeCost,lambda);
[minimum,va.selectedIndex]=min(va.cumulativeCost,[],1);
va.selectedIndex=va.selectedIndex(:);
va.selectedOmega=omega(va.selectedIndex);
va.selectedTime=time;
va.averageCost=bsxfun(@rdivide,va.cumulativeCost,1:nTime);
va.selectedAverageCost=minimum(:)./(1:nTime)';
va.q=-va.selectedAverageCost;
[~,va.rankOrder]=sortrows([-va.q,(1:nTime)'],[1 2]);
va.lambda=lambda; va.eta=0; va.epsilonP=epsilonP;
va.nonpositiveFraction=mean(wd(:)<=0);
va.clippedFraction=mean(normalized(:)<epsilonP);
va.tieCount=sum(bsxfun(@eq,va.cumulativeCost,minimum),1)';
va.tieColumns=sum(va.tieCount>1);
va.qualityTieCount=nTime-numel(unique(va.q));
va.qualityTieFraction=va.qualityTieCount/nTime;
end
