no_gran = {42
    37
    [46, 13]
    34
    24
    31
    37
    [13, 56];
    [34, 22];
    [19, 63];
    [23, 64];
    [39, 34];
    [27, 58]
    [23, 35]
    };

has_gran = {[24, 26]
    [56, 29]
    [19, 21, 26]
    [22, 17]
    [15, 23]
    [22, 17]
    [27, 14]
    [14, 17]
    [34, 22]
    [32, 38]
    [28, 26]
    [24, 36, 29]
    [29, 19]
    [28, 34]
    };

%% 
% 先确保计数变量存在
ng_cnt = cellfun(@numel, no_gran(:));
hg_cnt = cellfun(@numel, has_gran(:));

% swarm x
x1 = categorical(repmat("No granule", numel(ng_cnt), 1));
x2 = categorical(repmat("Has granule", numel(hg_cnt), 1));

figure('Color','w'); hold on
swarmchart(x1, ng_cnt, 36, 'filled', 'DisplayName','No granule');
swarmchart(x2, hg_cnt, 36, 'filled', 'DisplayName','Has granule');
ylabel('Divisions per cell'); box on

A = ng_cnt(:);  B = hg_cnt(:);
nA = numel(A);  nB = numel(B);

mA = mean(A); sA = std(A);
mB = mean(B); sB = std(B);

rng(1);
nboot = 10000;

boot_means_A = zeros(nboot,1);
for i = 1:nboot
    Ai = A(randi(nA, nA, 1));
    boot_means_A(i) = mean(Ai);
end
ciA = quantile(boot_means_A, [0.025 0.975]);

boot_means_B = zeros(nboot,1);
for i = 1:nboot
    Bi = B(randi(nB, nB, 1));
    boot_means_B(i) = mean(Bi);
end
ciB = quantile(boot_means_B, [0.025 0.975]);

loA = mA - ciA(1); hiA = ciA(2) - mA;
loB = mB - ciB(1); hiB = ciB(2) - mB;

errorbar(categorical("No granule"),  mA, loA, hiA, 'kd', ...
    'MarkerFaceColor','w', 'LineWidth',1.6, 'CapSize',10, ...
    'DisplayName','Mean ± 95% bootstrap CI');

errorbar(categorical("Has granule"), mB, loB, hiB, 'kd', ...
    'MarkerFaceColor','w', 'LineWidth',1.6, 'CapSize',10, ...
    'HandleVisibility','off');

legend({'No granule','Has granule','Mean ± 95% bootstrap CI (nboot = 10000)'}, 'Location','best');

ylim([0 4.5])
% 标题里说明：用的可视化与检验
title(sprintf('Divisions per cell by PHB granules (swarm + mean with 95%% CI; permutation test)'), 'Interpreter','tex', 'FontWeight','bold');


% 若你前面尚未把 cell 转为“每个细胞的分裂次数（1/2/3）”，先做这两行：
if ~exist('ng_cnt','var') || isempty(ng_cnt), ng_cnt = cellfun(@numel, no_gran(:)); end
if ~exist('hg_cnt','var') || isempty(hg_cnt), hg_cnt = cellfun(@numel, has_gran(:)); end
A = ng_cnt(:);  B = hg_cnt(:);

% 2×3 频数表（列：1/2/3 次）
edges = 0.5:1:3.5;                   % 把 1/2/3 落到各自的箱
cA = histcounts(A, edges);           % [n1, n2, n3] for No-granule
cB = histcounts(B, edges);           % [n1, n2, n3] for Has-granule
O  = [cA; cB];                       % 观测频数 2×3

% 卡方独立性检验
rowS = sum(O, 2);  colS = sum(O, 1);  N = sum(O, 'all');
E = rowS * colS / N;                  % 期望频数（独立假设下）
chi2 = sum((O - E).^2 ./ max(E, eps), 'all');
df = (size(O,1)-1) * (size(O,2)-1);
p_chi = 1 - chi2cdf(chi2, df);

% （可选）合并成 2×2 做 Fisher：1 次 vs ≥2 次
T2 = [sum(A==1), sum(A>=2); sum(B==1), sum(B>=2)];
p_fisher = NaN;
% try
%     [~, p_fisher] = fishertest(T2);  % 需要统计工具箱；没有则保持 NaN
% catch
%     % 若无工具箱，可忽略或自己实现精确检验
% end

% ——把结果写进图里（和你之前风格一致的白底文字框）——
ax = gca;
txt_chi = sprintf(['$\\chi^2(%d)=%.2f,\\; p=%.3g$,\\ ', ...
                   'Counts (1/2/3): No %d/%d/%d; Has %d/%d/%d'], ...
                   df, chi2, p_chi, cA(1), cA(2), cA(3), cB(1), cB(2), cB(3));

if ~isnan(p_fisher)
    txt_chi = sprintf('%s \\\\ Fisher (1 vs $\\geq$2): $p=%.3g$', txt_chi, p_fisher);
end

text(ax, 0.92, 0.12, txt_chi, ...
    'Units','normalized', ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','top', ...
    'BackgroundColor','w', ...
    'EdgeColor',[0.2 0.2 0.2], ...
    'Margin', 6, ...
    'Interpreter','latex', ...
    'FontWeight','bold');




%% 

% —— 把每个细胞的时间拉直并合并 —— 
no_gran  = cellfun(@(x) sort(x(:)), no_gran,  'UniformOutput', false);
has_gran = cellfun(@(x) sort(x(:)), has_gran, 'UniformOutput', false);

times_no  = vertcat(no_gran{:});   % No granule 组所有分裂时间
times_has = vertcat(has_gran{:});  % Has granule 组所有分裂时间

% —— 自动选择共同的 bin（用合并后的数据来定），避免两组各自的 bin 不一致 —— 
all_times = [times_no; times_has];
if isempty(all_times)
    error('两组时间数据都是空的。');
end
% Freedman–Diaconis 规则选 bin 宽，兜底防止过窄/过宽
bw_fd = 2*iqr(all_times)/max(1,numel(all_times)^(1/3));
if ~isfinite(bw_fd) || bw_fd<=0
    bw_fd = max(1, std(all_times)/5);  % 兜底
end
lo = min(all_times); hi = max(all_times);
edges = lo:bw_fd:(hi+bw_fd);
if numel(edges) < 3
    edges = linspace(lo, hi, 20);
end

% —— 画重叠直方图（概率密度归一化 + 半透明）——
figure('Name','All division times (pooled)','Color','w');
hold on
h1 = histogram(times_no,  edges, 'Normalization','pdf', 'DisplayName','No granule');
h2 = histogram(times_has, edges, 'Normalization','pdf', 'DisplayName','Has granule');
% 半透明方便叠加对比
if isprop(h1,'FaceAlpha')
    h1.FaceAlpha = 0.5;
    h2.FaceAlpha = 0.5;
end
xlabel('分裂时间'); ylabel('概率密度');
title('两组所有分裂时间的分布（合并后叠加）');
legend('Location','best'); box on



%% 
% ==== Figure 2: All division times (boxplot) ====
data  = [times_no; times_has];
group = [repmat({'No granule'}, numel(times_no), 1);
         repmat({'Has granule'}, numel(times_has), 1)];

figure('Name','Division times (boxplot)','Color','w');
boxplot(data, group, 'Notch','on', 'Whisker',1.5);
ylabel('Division time');
title('Boxplot of all division times (pooled)');
% === Welch's t-test (unequal variances) & annotate ===
if exist('ttest2','file')~=2
    warning('ttest2 is not available. Need Statistics and Machine Learning Toolbox.');
else
    % Welch t-test on pooled division times
    [~, p, ci, stats] = ttest2(times_no, times_has, 'Vartype','unequal');  % ci for mean difference (no - has)
    dmu = mean(times_no) - mean(times_has);

    % Put summary text on the existing boxplot
    ax = gca;
    yl = ylim(ax);
    rngY = yl(2) - yl(1);
    ytxt = yl(2) + 0.12*rngY;   % text y-position (a bit above top)
    ybar = yl(2) + 0.06*rngY;   % bar y-position

    % Make room above for text/bar
    ylim([yl(1), yl(2) + 0.2*rngY]);

    % significance stars
    if p < 1e-3
        stars = '***';
    elseif p < 1e-2
        stars = '**';
    elseif p < 0.05
        stars = '*';
    else
        stars = 'n.s.';
    end

    % Draw a significance bar between group 1 and 2
    hold on
    % plot([1 1 2 2], [ybar-0.005*rngY ybar ybar ybar-0.005*rngY], 'k-', 'LineWidth', 1.2)
    % text(1.5, ybar + 0.01*rngY, stars, 'HorizontalAlignment','center', 'FontWeight','bold')

    % Summary text box
    txt = sprintf(['Welch t-test: t(%0.1f) = %0.3f, p = %0.3g\n' ...
                   '95%% CI for \\Delta\\mu (no - has): [%0.3g, %0.3g]\n' ...
                   '\\Delta\\mu = %0.3g'], stats.df, stats.tstat, p, ci(1), ci(2), dmu);
    text(1.7, ytxt, txt, ...
        'HorizontalAlignment','center', 'VerticalAlignment','top', ...
        'BackgroundColor','w', 'Margin',6, 'EdgeColor',[0.8 0.8 0.8], 'FontName','Arial');
end

grid on

