function verify_wdr_va
% 数值验证：在 MATLAB 中点击运行，不保存任何输出文件。
scriptDir=fileparts(mfilename('fullpath')); addpath(fullfile(fileparts(scriptDir),'modules'));
u=[1 4 4;4 1 4]; [d,pred]=wdr_va_recurrence(u,1);
assert(isequal(d,[1 5 8;4 3 7])); assert(isequal(pred(:,2),[1;1]));
[v,k]=min(d,[],1); avg=v./(1:3); q=-avg;
[~,order]=sortrows([-q(:),(1:3)'],[1 2]);
assert(isequal(k,[1 2 2]) && max(abs(avg-[1 1.5 7/3]))<1e-12);
assert(isequal(order,[1;2;3]));

% 对小矩阵穷举每条路径，比较每个时刻和终点；也检查 lambda=0。
u=[1 2 4 0;2 0 3 1;3 1 0 2];
for lambda=[0 0.7]
    d=wdr_va_recurrence(u,lambda); nf=size(u,1);
    for n=1:size(u,2)
        brute=inf(nf,1);
        for code=0:nf^n-1
            value=code; path=zeros(1,n);
            for j=1:n, path(j)=mod(value,nf)+1; value=floor(value/nf); end
            linear=sub2ind(size(u),path,1:n);
            candidate=sum(u(linear))+lambda*sum(diff(path).^2);
            brute(path(end))=min(brute(path(end)),candidate);
        end
        assert(max(abs(d(:,n)-brute))<1e-10);
    end
end

% 全局正值归一化、零值和负值截断。
va=wdr_va_select([1 2;0 -1],[0;1],[0;1],0,1e-12);
expected=-log([1/3 2/3;1e-12 1e-12]);
assert(max(abs(va.nodeCost(:)-expected(:)))<1e-12);
assert(va.nonpositiveFraction==0.5 && va.clippedFraction==0.5);
failed=false; try, wdr_va_select(-ones(2,2),[0;1],[0;1],0,1e-12); catch ex, failed=strcmp(ex.identifier,'wdr:VaNonpositive'); end; assert(failed);
failed=false; try, wdr_va_select([1 NaN;1 1],[0;1],[0;1],0,1e-12); catch, failed=true; end; assert(failed);
failed=false; try, wdr_va_select(ones(3,2),[0;1;3],[0;1],0,1e-12); catch ex, failed=strcmp(ex.identifier,'wdr:VaFrequency'); end; assert(failed);
failed=false; try, wdr_va_select(ones(3,2),[0;1;2],[0;0],0,1e-12); catch ex, failed=strcmp(ex.identifier,'wdr:VaTime'); end; assert(failed);
va=wdr_va_select(ones(3,3),(1:3)',(1:3)',0,1e-12);
assert(all(va.selectedIndex==1) && va.tieColumns==3 && isequal(va.rankOrder,(1:3)'));
va=wdr_va_select([10 1;1 1;1 10],(1:3)',(1:2)',100,1e-12);
assert(va.selectedIndex(2)==1); % 第二列选点不同于 WD 最大峰第 3 行。

cfg=wdr_config('quick'); rng(42,'twister'); signal=wdr_generate_lfm(cfg,-5,[220 -70]);
[wd,omega,t]=wdr_distribution(signal.x,cfg); ridge=wdr_ridge(wd,omega);
rng(3,'twister'); [est,va,schedule]=wdr_va_prosac_estimate(signal.x,cfg.t,t,wd,omega,[5 10],cfg,0.032,1e-12,10);
assert(isequal(sort(va.rankOrder),(1:numel(t))'));
assert(all(schedule.pairs(1,:)~=schedule.pairs(2,:)));
assert(isequal(size(est.omega),[numel(t),2]));
for j=1:10
    i1=schedule.pairs(1,j); i2=schedule.pairs(2,j);
    a=(va.selectedOmega(i2)-va.selectedOmega(i1))/(t(i2)-t(i1));
    b=va.selectedOmega(i1)-a*t(i1);
    manual=abs(sum(signal.x.*exp(-1j*(a*cfg.t.^2/2+b*cfg.t))));
    assert(abs(manual-est.candidateScores(j))<1e-8);
end
assert(isequal(ridge,wdr_ridge(wd,omega)));
cfg.trials=2; cfg.snrDb=[-5 5]; cfg.mList=[5 10]; cfg.showProgress=false;
old=wdr_compare_prosac(cfg,10); compared=wdr_compare_va(cfg,0.032,1e-12,10);
assert(isequal(old.gaussianError,compared.baselineError));
assert(all(isfinite(compared.vaProsacMse(:))));
assert(all(compared.vaProsacTotalSeconds(:)>=compared.vaSeconds(:)));
fprintf('通过：手算、穷举、边界、选点、抽样、相干评分及原基线对照。\n');
end
