A = total_cell_data_collect;   % 行=特征，列=样本/周期

% 差值：第12行 - 第11行
d = (A(12,:) - A(11,:))/6;

% 分组：按第3行是否为0
mask_zero3    = (A(3,:) == 0);
mask_nonzero3 = ~mask_zero3;

% 过滤无效值
valid = isfinite(d);
g1 = d(valid & mask_zero3);      % 第3行为0的差值
g2 = d(valid & mask_nonzero3);   % 第3行非零的差值

%% 
% 确保 g1, g2 已由你上一步构造好；这里再做一次健壮性过滤
g1 = g1(isfinite(g1));
g2 = g2(isfinite(g2));

% 均值差（按 MATLAB ttest2 的定义，差=mean(g1)-mean(g2)）
delta = mean(g1) - mean(g2);

% Welch 双样本 t 检验（95% CI）
[H,P,CI,STATS] = ttest2(g1, g2, 'Vartype','unequal', 'Alpha',0.05);

% 输出结果
fprintf('Δμ = %.4f (h), 95%% CI [%.4f, %.4f]\n', delta, CI(1), CI(2));
fprintf('Welch t-test: t(%0.1f) = %.3f, p = %.3g, H = %d\n', STATS.df, STATS.tstat, P, H);


%% 

% 统一的 bin 边界（示例：固定 bin 宽度）
binWidth = 0.2;  % 想要别的宽度就改这里
allvals = [g1(:); g2(:)];
lo = floor(min(allvals)/binWidth)*binWidth;
hi = ceil(max(allvals)/binWidth)*binWidth;
edges = lo:binWidth:hi;


% 绘图：两组直方图叠在同一张图上
figure; hold on;
h1 = histogram(g1, 'BinEdges', edges, 'Normalization', 'probability');
h2 = histogram(g2, 'BinEdges', edges, 'Normalization', 'probability');

% 让两组直方图更易区分：透明度 & 取消边框
h1.FaceAlpha = 0.5;  h2.FaceAlpha = 0.5;
% h1.EdgeColor = 'none'; h2.EdgeColor = 'none';
xlim([0 10])
xlabel('Row12 - Row11');
ylabel('Fraction of columns');  % 归一化为概率后的纵轴
title('Overlayed histograms of doubling time grouped by granule existance');
legend([h1 h2], {'granule-free daughters', 'granule-bearing daughters'}, 'Location', 'best');
% 可选：在当前图上加注释（或你新建一幅图后再标）
txt = sprintf('\\Delta\\mu=%.4f, 95%% CI [%.4f, %.4f]; t(%.1f)=%.3f, p=%.3g', ...
              delta, CI(1), CI(2), STATS.df, STATS.tstat, P);
annotation('textbox',[0.35 0.52 0.5 0.08],'String',txt,'EdgeColor','none','HorizontalAlignment','center');
grid on; hold off;


%% 
A = total_cell_data_collect;     % 行=特征，列=样本/周期

% 数据
x = A(3,:);                                % 横轴：第3行
y = (A(12,:) - A(11,:)) / 6;               % 纵轴：差值/6

% 清洗无效值
mask = isfinite(x) & isfinite(y);
x = x(mask);
y = y(mask);

% 画散点
figure; hold on;
sc = scatter(x, y, 8, 'filled');    % 点小一些
grid on;
xlabel('daughter granule size');
ylabel('doubling time (h)');
title('Granule size and division rate with linear fit');

% —— 线性拟合（优先用 fitlm 拿到 slope/CI/R^2/p）——
haveStatsTBX = exist('fitlm','file') == 2;
if haveStatsTBX
    mdl = fitlm(x', y');                  % 注意: 列向量
    m   = mdl.Coefficients.Estimate(2);   % slope
    b   = mdl.Coefficients.Estimate(1);   % intercept
    CI  = coefCI(mdl, 0.05);              % 95% CI
    slope_CI = CI(2,:);                   % slope 的CI
    R2  = mdl.Rsquared.Ordinary;
    pval = mdl.Coefficients.pValue(2);
else
    % 备用：没有统计工具箱时，用 polyfit
    p = polyfit(x, y, 1);
    m = p(1); b = p(2);
    % 手算 R^2
    yhat = polyval(p, x);
    SSres = sum((y - yhat).^2);
    SStot = sum((y - mean(y)).^2);
    R2 = 1 - SSres/SStot;
    pval = NaN; slope_CI = [NaN NaN];
end

% 画回归直线
xx = linspace(min(x), max(x), 200);
yy = m*xx + b;
ln = plot(xx, yy, 'LineWidth', 2);

% 在图上标注 slope（以及 R^2 / p / CI）
if haveStatsTBX
    txt = sprintf('slope = %.4f (95%% CI [%.4f, %.4f])\nR^2 = %.3f, p = %.3g', ...
                  m, slope_CI(1), slope_CI(2), R2, pval);
else
    txt = sprintf('slope = %.4f\\nR^2 = %.3f', m, R2);
end
legend([sc ln], {'data','linear fit'}, 'Location','best');
annotation('textbox',[0.15 0.62 0.7 0.12], 'String', txt, ...
           'EdgeColor','none', 'HorizontalAlignment','center');

hold off;
