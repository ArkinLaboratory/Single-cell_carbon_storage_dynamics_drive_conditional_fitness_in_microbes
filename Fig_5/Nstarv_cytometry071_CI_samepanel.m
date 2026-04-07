%% 原始数据（同之前）
% --- 原始数据：WT 三个样本 ---
time = ["12:24";
        "12:38";
        "13:55";
        "14:49";
        "15:44";
        "16:48";
        "17:48";
        "18:50";
        "19:50";
        "20:50";
        "21:50"];

WT1 = [15538;17317;15724;19145;18819;19032;23588;25301;29812;36298;43076];
WT2 = [27976;27804;30468;30236;34583;34330;46150;46179;48506;70770;90228];
WT3 = [28055;25252;26494;30320;30947;35427;45633;48841;56496;70225;95094];

% --- 原始数据：ΔPhaC 三个样本 ---
dPha1 = [17838;17446;16046;18188;15702;16889;20571;19571;20282;25021;30087];
dPha2 = [17795;17612;17463;17695;18837;17691;22404;21841;20068;26492;32671];
dPha3 = [11958;11060;11016;12151;12782;12212;13892;14663;13036;14997;25791];

% --- 参考组(ref) ---
time_ref = ["11:07";"12:03";"13:07";"14:02";"15:08";"16:08";"17:06";"18:05";"19:05"];
WT_ref  = [72502;73592;71347;75584;77748;80931;66391;84840;86688];
dPha_ref= [92434;87593;75406;65759;69920;69560;77082;79288;84832];

WT1_norm   = WT1   ./ WT1(1);     dPha1_norm = dPha1 ./ dPha1(1);
WT2_norm   = WT2   ./ WT2(1);     dPha2_norm = dPha2 ./ dPha2(1);
WT3_norm   = WT3   ./ WT3(1);     dPha3_norm = dPha3 ./ dPha3(1);
WT_ref_norm   = WT_ref   ./ WT_ref(1);
dPha_ref_norm = dPha_ref ./ dPha_ref(1);

normM1    = [WT1_norm, dPha1_norm];
normM2    = [WT2_norm, dPha2_norm];
normM3    = [WT3_norm, dPha3_norm];
normM_ref = [WT_ref_norm, dPha_ref_norm];

res1 = [cellstr(time),     num2cell(normM1)];
res2 = [cellstr(time),     num2cell(normM2)];
res3 = [cellstr(time),     num2cell(normM3)];
res_ref = [cellstr(time_ref), num2cell(normM_ref)];

%% 清理环境
clear; clc;

%% 1. 原始数据（归一化后的）
t1 = [0.000, 0.233, 1.517, 2.417, 3.333, 4.400, 5.400, 6.433, 7.433, 8.433, 9.433];
t2 = t1;
t3 = t1;
t4 = [0.000, 0.933, 2.000, 2.917, 4.017, 5.017, 5.983, 6.967, 7.967];

WT1  = [1.000, 1.114, 1.012, 1.232, 1.211, 1.225, 1.518, 1.628, 1.919, 2.336, 2.772];
WT2  = [1.000, 0.994, 1.089, 1.081, 1.236, 1.227, 1.650, 1.651, 1.734, 2.530, 3.225];
WT3  = [1.000, 0.900, 0.944, 1.081, 1.103, 1.263, 1.627, 1.741, 2.014, 2.503, 3.390];
WT4  = [1.000, 1.015, 0.984, 1.043, 1.072, 1.116, 0.916, 1.170, 1.196];

DP1  = [1.000, 0.978, 0.900, 1.020, 0.880, 0.947, 1.153, 1.097, 1.137, 1.403, 1.687];
DP2  = [1.000, 0.990, 0.981, 0.994, 1.059, 0.994, 1.259, 1.227, 1.128, 1.489, 1.836];
DP3  = [1.000, 0.925, 0.921, 1.016, 1.069, 1.021, 1.162, 1.226, 1.090, 1.254, 2.157];
DP4  = [1.000, 0.948, 0.816, 0.711, 0.756, 0.753, 0.834, 0.858, 0.918];

%% 2. 合并 1–3 组，用于拟合
t_exp  = repmat(t1', 3, 1);            % 33×1
WT_exp = [WT1'; WT2'; WT3'];           % 33×1
DP_exp = [DP1'; DP2'; DP3'];           % 33×1
N      = numel(t_exp);

%% 3. 带 lag 的分段平滑指数模型（函数形式，返回 [WT; DP] 列向量）
smoothstep = @(x,lambda,b) 0.5 + 0.5*tanh((x - lambda)/b);
fixed_b    = 2.45;

p0   = [max(WT_exp), 0.2, 4,  max(DP_exp), 0.2, 4];
lb   = [1, 0.01, 0,   1,        0.01, 0];
ub   = [5, 1, max(t1), 3,       1,    max(t1)];
opts = optimoptions('lsqcurvefit','MaxIter',1000,'TolFun',1e-10,'Display','off');

[params_lag,~,res_lag,~,~,~,J_lag] = lsqcurvefit(@modlag_vec, ...
    p0, t_exp, [WT_exp; DP_exp], lb, ub, opts);

%% 4. 不带 lag 的纯指数模型（函数形式，返回 [WT; DP] 列向量）
p0_nl = [1,0.2,  1,0.2];
lb_nl = [0,0,    0,0];
ub_nl = [5,1,    5,1];

[params_nl,~,res_nl,~,~,~,J_nl] = lsqcurvefit(@modnl_vec, ...
    p0_nl, t_exp, [WT_exp; DP_exp], lb_nl, ub_nl, opts);

%% 5. 参考组：线性模型（便于获得预测区间）
mdlRef_WT = fitlm(t4', WT4');   % linearmodel + predict → y & CI
mdlRef_DP = fitlm(t4', DP4');

%% 6. 生成拟合曲线与 95% 置信带
alpha = 0.05;
xFit  = linspace(min(t_exp), max(t_exp), 300).';   % 列向量
nF    = numel(xFit);

% lag CI
[yLag_all, dLag_all] = nlpredci(@modlag_vec, xFit, params_lag, res_lag, ...
                                'Jacobian', J_lag, 'Alpha', alpha);
yLag_WT = yLag_all(1:nF);      dLag_WT = dLag_all(1:nF);
yLag_DP = yLag_all(nF+1:end);  dLag_DP = dLag_all(nF+1:end);
loLag_WT = yLag_WT - dLag_WT;  hiLag_WT = yLag_WT + dLag_WT;
loLag_DP = yLag_DP - dLag_DP;  hiLag_DP = yLag_DP + dLag_DP;

% no-lag CI
[yNL_all, dNL_all] = nlpredci(@modnl_vec, xFit, params_nl, res_nl, ...
                              'Jacobian', J_nl, 'Alpha', alpha);
yNL_WT = yNL_all(1:nF);        dNL_WT = dNL_all(1:nF);
yNL_DP = yNL_all(nF+1:end);    dNL_DP = dNL_all(nF+1:end);
loNL_WT = yNL_WT - dNL_WT;     hiNL_WT = yNL_WT + dNL_WT;
loNL_DP = yNL_DP - dNL_DP;     hiNL_DP = yNL_DP + dNL_DP;

% ref 线性预测 + CI
xRef = linspace(min(t4), max(t4), 300).';
[yRef_WT, yci_WT] = predict(mdlRef_WT, xRef);   % 返回 y 与 CI（可指定 Alpha/Prediction）.
[yRef_DP, yci_DP] = predict(mdlRef_DP, xRef);   % :contentReference[oaicite:1]{index=1}

%% 7. 在训练点上做预测并计算性能指标（保留你的统计输出）
Ylag_all     = modlag_vec(params_lag, t_exp);
predLagWT    = Ylag_all(1:N);
predLagDP    = Ylag_all(N+1:end);

Ynl_all      = modnl_vec(params_nl, t_exp);
predNoLagWT  = Ynl_all(1:N);
predNoLagDP  = Ynl_all(N+1:end);

SSres_lag_WT = sum((WT_exp - predLagWT).^2);
SStot_WT     = sum((WT_exp - mean(WT_exp)).^2);
R2_lag_WT    = 1 - SSres_lag_WT/SStot_WT;

SSres_lag_DP = sum((DP_exp - predLagDP).^2);
SStot_DP     = sum((DP_exp - mean(DP_exp)).^2);
R2_lag_DP    = 1 - SSres_lag_DP/SStot_DP;

SSres_nl_WT  = sum((WT_exp - predNoLagWT).^2);
R2_nl_WT     = 1 - SSres_nl_WT/SStot_WT;
SSres_nl_DP  = sum((DP_exp - predNoLagDP).^2);
R2_nl_DP     = 1 - SSres_nl_DP/SStot_DP;

R2c_lag_WT = corr(WT_exp, predLagWT)^2;
R2c_lag_DP = corr(DP_exp, predLagDP)^2;
R2c_nl_WT  = corr(WT_exp, predNoLagWT)^2;
R2c_nl_DP  = corr(DP_exp, predNoLagDP)^2;

fprintf('\n==== Model Performance ====\n');
fprintf('WT (lag)     SSres=%.3f  SStot=%.3f  R2=%.3f  Corr2=%.3f\n', SSres_lag_WT, SStot_WT, R2_lag_WT, R2c_lag_WT);
fprintf('WT (no-lag)  SSres=%.3f            R2=%.3f  Corr2=%.3f\n',   SSres_nl_WT,            R2_nl_WT,  R2c_nl_WT);
fprintf('ΔP (lag)     SSres=%.3f  SStot=%.3f  R2=%.3f  Corr2=%.3f\n', SSres_lag_DP, SStot_DP, R2_lag_DP, R2c_lag_DP);
fprintf('ΔP (no-lag)  SSres=%.3f            R2=%.3f  Corr2=%.3f\n\n', SSres_nl_DP,            R2_nl_DP,  R2c_nl_DP);

%% 8. 单面板绘图（OpenGL 以支持透明度）
% —— 矢量可编辑版（无透明）——
% === 矢量可编辑 + 柔和“假透明”CI（无FaceAlpha） ===
figure('Renderer','painters'); hold on;

WTColor = [0 0.4470 0.7410];    % 蓝
DPColor = [0.8500 0.3250 0.0980]; % 橙

% 辅助函数：把颜色朝白色拉 f∈(0,1)，f越大越浅
lighten = @(c,f) c*(1-f) + f;   % e.g., lighten([r g b], 0.88)

WTshade = lighten(WTColor, 0.90);   % 更浅，更像半透明
DPshade = lighten(DPColor, 0.90);

% —— 保证 x 单调唯一，避免自交（CI带边缘更顺）——
[xFill, iuniq] = unique(xFit(:).');        % 行向量、唯一
loWT = loLag_WT(iuniq).';   hiWT = hiLag_WT(iuniq).';
loDP = loLag_DP(iuniq).';   hiDP = hiLag_DP(iuniq).';

[xRefFill, iri] = unique(xRef(:).');
loWTref = yci_WT(iri,1).';  hiWTref = yci_WT(iri,2).';
loDPref = yci_DP(iri,1).';  hiDPref = yci_DP(iri,2).';

% —— 数据点（同组形状对应）——
scatter(t1,WT1,55,WTColor,'filled','o','DisplayName','WT G1');
scatter(t2,WT2,55,WTColor,'filled','s','DisplayName','WT G2');
scatter(t3,WT3,55,WTColor,'filled','^','DisplayName','WT G3');
scatter(t4,WT4,55,WTColor,'d','MarkerFaceColor','w','LineWidth',1.2,'DisplayName','WT ref');

scatter(t1,DP1,50,DPColor,'filled','o','DisplayName','\DeltaPhaC G1');
scatter(t2,DP2,50,DPColor,'filled','s','DisplayName','\DeltaPhaC G2');
scatter(t3,DP3,50,DPColor,'filled','^','DisplayName','\DeltaPhaC G3');
scatter(t4,DP4,50,DPColor,'d','MarkerFaceColor','w','LineWidth',1.2,'DisplayName','\DeltaPhaC ref');

% —— CI“假透明”阴影（不透明浅色 patch，矢量安全）——
pWT  = fill([xFill fliplr(xFill)], [loWT fliplr(hiWT)], WTshade, 'EdgeColor','none', 'DisplayName','WT 95% CI');
pDP  = fill([xFill fliplr(xFill)], [loDP fliplr(hiDP)], DPshade, 'EdgeColor','none', 'DisplayName','\DeltaPhaC 95% CI');
pWT2 = fill([xRefFill fliplr(xRefFill)], [loWTref fliplr(hiWTref)], WTshade, 'EdgeColor','none','HandleVisibility','off');
pDP2 = fill([xRefFill fliplr(xRefFill)], [loDPref fliplr(hiDPref)], DPshade, 'EdgeColor','none','HandleVisibility','off');
uistack([pWT pDP pWT2 pDP2],'bottom');   % 阴影放底层

% —— 拟合曲线（设置端点圆角可读性更好；AI 中也能改）——
plot(xFit, yLag_WT, '-', 'Color', WTColor, 'LineWidth', 2.2, 'DisplayName','WT fit (lag)');
plot(xFit, yNL_WT,  ':', 'Color', WTColor, 'LineWidth', 2.0, 'DisplayName','WT fit (no-lag)');
plot(xFit, yLag_DP, '-', 'Color', DPColor,'LineWidth', 2.2, 'DisplayName','\DeltaPhaC fit (lag)');
plot(xFit, yNL_DP,  ':', 'Color', DPColor,'LineWidth', 2.0, 'DisplayName','\DeltaPhaC fit (no-lag)');

% —— ref 线性拟合 —— 
plot(xRef, yRef_WT, '--', 'Color', WTColor, 'LineWidth',2.0, 'DisplayName','WT lin ref');
plot(xRef, yRef_DP, '--', 'Color', DPColor,'LineWidth',2.0, 'DisplayName','\DeltaPhaC lin ref');

xlabel('Time (h)'); ylabel('Normalized value');
title('WT & \DeltaPhaC vs Time (single panel, vector-safe CI)');
ylim([0.5 3.5]); grid on;
legend('Location','bestoutside');

% —— 导出为矢量（嵌入字体；AI 可逐元素编辑）——
exportgraphics(gcf,'figure_vector.pdf','ContentType','vector');   % 或 print -dpdf -painters

hold off;


%% —— 本脚本的局部函数 —— 
% 说明：MATLAB 支持在脚本中定义局部函数（R2024a 以前需放在文件结尾）。:contentReference[oaicite:2]{index=2}

function y = modlag_vec(p, x)
% 带 lag 的分段平滑指数模型；返回 [WT(:); DP(:)] 列向量
x = x(:);
fixed_b   = 2.45;
smoothstep = @(z,lambda,b) 0.5 + 0.5*tanh((z - lambda)/b);
WT = 1 + (p(1).*exp(p(2).*(x - p(3))) - 1) .* smoothstep(x, p(3), fixed_b);
DP = 1 + (p(4).*exp(p(5).*(x - p(6))) - 1) .* smoothstep(x, p(6), fixed_b);
y  = [WT; DP];
end

function y = modnl_vec(p, x)
% 不带 lag 的纯指数模型；返回 [WT(:); DP(:)] 列向量
x  = x(:);
WT = p(1).*exp(p(2).*x);
DP = p(3).*exp(p(4).*x);
y  = [WT; DP];
end
