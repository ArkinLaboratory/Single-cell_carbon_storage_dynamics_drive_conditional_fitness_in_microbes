%% 0) 统计每个 bin 的死亡数与样本数（沿用你的写法）
step = 0.05;
frac1 = round(total_cell_data_collect(16,:) / step) * step;
total_length = numel(frac1);

frac = [frac1; zeros(1,total_length)];
for i = 1:total_length
    if i==total_length || total_cell_data_collect(10,i) ~= total_cell_data_collect(10,i+1)
        frac(2,i) = frac(2,i) + 1;     % 你定义的“死亡”事件记+1
    end
end

ux = unique(frac(1,:));                  % 所有分箱
uniq_frac = zeros(5,numel(ux));          % [x; dead; total; rate; se]
uniq_frac(1,:) = ux;
for k = 1:numel(ux)
    cols = (frac(1,:)==ux(k));
    dead_cnt  = sum(frac(2,cols));
    total_cnt = sum(cols);
    p = dead_cnt/total_cnt;
    uniq_frac(2,k)=dead_cnt; uniq_frac(3,k)=total_cnt;
    uniq_frac(4,k)=p;        uniq_frac(5,k)=sqrt(p*(1-p)/max(total_cnt,1));
end

%% === 4参数Logistic（4PL）拟合：仅用蓝点 ===
x_pts = uniq_frac(1, :)';   % 蓝点 x
y_pts = uniq_frac(4, :)';   % 蓝点 y（死亡率）

% 1) 基本清洗与排序
mask = isfinite(x_pts) & isfinite(y_pts);
x = x_pts(mask);  y = y_pts(mask);
[x, idx] = sort(x);  y = y(idx);

% 2) 4PL 模型与残差（等权；加一轮 Huber 稳健权重可选）
log4 = @(p,xx) p(1) + (p(2)-p(1)) ./ (1 + exp(-p(3).*(xx - p(4)))); % p=[L,U,k,x0]
% 初值：用数据的下/上分位做平台，x0 先放在0.7附近（或x的中位），k给个中等斜率
L0  = max(0,  min(y));
U0  = min(1,  max(y));
x0_ = median(x);                 % 如果你期望0.7，可直接设 x0_=0.7
x0_ = max(min(x0_, max(x)-0.05*(max(x)-min(x))), min(x)+0.05*(max(x)-min(x)));
k0  = 20;                        % 陡峭度初值（10–40 常见）
p0  = [L0, U0, k0, x0_];

% 边界：保证 0<=L<U<=1，k>=0，x0 在范围内
lb = [0,   0.3,  0,   min(x)];
ub = [0.7, 1.0,  200, max(x)];

% 3) 先等权拟合，再做一轮稳健权重
opts = optimoptions('lsqnonlin','Display','off');

% 等权一轮
resfun = @(p) log4(p,x) - y;
p_hat  = lsqnonlin(resfun, p0, lb, ub, opts);

% Huber 稳健一轮（可选，但建议保留）
r  = log4(p_hat, x) - y;
c  = 1.345 * median(abs(r - median(r)) + eps);
w  = ones(size(r)); a = abs(r)/c; w(a>1) = 1./a(a>1);
resfun_w = @(p) w .* (log4(p,x) - y);
p_hat = lsqnonlin(resfun_w, p_hat, lb, ub, opts);

% 4) 预测与绘图
xg = linspace(min(x), max(x), 400)';
yg = log4(p_hat, xg);
L=p_hat(1); U=p_hat(2); k=p_hat(3); x0=p_hat(4);

figure; hold on; box on
% 原蓝点（含误差棒）
errorbar(uniq_frac(1,:), uniq_frac(4,:), uniq_frac(5,:), 'o', ...
    'MarkerSize', 6, 'MarkerEdgeColor','blue', 'MarkerFaceColor','blue', ...
    'CapSize', 10, 'DisplayName','Binned rate');

% 4PL 拟合曲线
plot(xg, yg, 'LineWidth', 2.2, 'Color', [0.2 0.6 0.9], ...
    'DisplayName','4-parameter Logistic fit');

s_low  = (3 - sqrt(3))/6;    % ~0.2113 -> y'' 最大（正）
s_high = (3 + sqrt(3))/6;    % ~0.7887 -> y'' 最小（负）

x_ddmax  = x0 - (1/k) * log((1 - s_low)  / s_low);   % 斜率增加最快（左侧）
x_ddmin  = x0 - (1/k) * log((1 - s_high) / s_high);  % 斜率减少最快（右侧）

% 画竖线（可选：只画左侧那个“增加最快”的）
plot([x_ddmax x_ddmax], [0 1], 'k:', 'LineWidth', 1.4, 'DisplayName','max d^2y/dx^2');
text(x_ddmax, 0.92, sprintf('  d^2 max = %.3f', x_ddmax), ...
     'VerticalAlignment','bottom', 'FontAngle','italic');

% 拐点（斜率最大处，恰好是 x0）竖线
ylim([0 1]); xlim([min(xg) max(xg)]);


% plot([x0 x0], [0 1], 'k--', 'LineWidth', 1.2, 'DisplayName','Inflection x_0');
% text(x0, 1, sprintf('  x_0 = %.3f', x0), 'VerticalAlignment','top');
% 
% % 注：斜率最大值 = (U-L)*k/4，可一起标注
% smax = (U - L) * k / 4;
% text(x0, 0.05, sprintf('  slope_{max}=%.2f', smax), 'VerticalAlignment','bottom');

xlabel('PHB fraction'); ylabel('Death rate');
title('4-parameter Logistic fit (points-only)');
legend('Location','best');


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
xlabel('Row 3');
ylabel('(Row12 - Row11) / 6');
title('Scatter with linear fit');

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
    txt = sprintf('slope = %.4f (95%% CI [%.4f, %.4f])\\nR^2 = %.3f, p = %.3g', ...
                  m, slope_CI(1), slope_CI(2), R2, pval);
else
    txt = sprintf('slope = %.4f\\nR^2 = %.3f', m, R2);
end
legend([sc ln], {'data','linear fit'}, 'Location','best');
annotation('textbox',[0.15 0.82 0.7 0.12], 'String', txt, ...
           'EdgeColor','none', 'HorizontalAlignment','center');

hold off;
