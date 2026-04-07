%% Plot granule fraction distribution (log-y) with minimal-model overlay
% Assumes:
%   - total_cell_data_collect exists in workspace (matrix)
%   - row 16 contains granule fraction values used in the original plot
%
% Output:
%   - Histogram of granule fraction (log y)
%   - Overlaid shifted-exponential (optionally truncated to [0,1]) curve
%
% Notes:
%   - "shifted exponential" is the distribution predicted by:
%       exponential age distribution + linear growth map (minimal model logic)
%   - We estimate x0 (shift) from a low percentile to be robust to noise.

clearvars -except total_cell_data_collect; 

%% ---------- 1) Extract data ----------
if ~exist('total_cell_data_collect','var')
    error('total_cell_data_collect not found in workspace. Load your .mat file first.');
end

% Try to interpret "row 16" robustly:
[m,n] = size(total_cell_data_collect);
if m >= 16
    x = total_cell_data_collect(16,:);
elseif n >= 16
    x = total_cell_data_collect(:,16);
else
    error('Neither dimension has length >= 16. Please verify where the granule fraction vector is stored.');
end
x = x(:);

% Clean
x = x(isfinite(x));

% If your statement is "within granule-bearing cells", filter >0
x = x(x >= 0);

% If you know fraction is bounded [0,1], optionally enforce it (comment out if not)
x = x(x <= 1);

if numel(x) < 50
    warning('Very few data points after filtering (%d). Check filters (>0, <=1).', numel(x));
end

%% ---------- 2) Choose histogram settings ----------
nbins = 80;                 % adjust to match original plot
use_pdf = true;             % true: normalize to PDF; false: raw counts
log_y = true;

% For overlay curve, we will plot on the same scale as histogram
% If PDF: overlay PDF; if counts: overlay scaled counts

%% ---------- 3) Fit shifted-exponential model ----------
% Choose shift x0 as a low percentile (robust against tiny noisy values)
x0 = prctile(x, 5);         % you can try 1 or 10 depending on your data
x_fit = x(x >= x0);

% MLE for exponential on (x - x0): lambda = 1/mean(x-x0)
dx = x_fit - x0;
dx = dx(dx > 0);
lambda = 1 / mean(dx);

% If fraction is bounded above by 1, you may prefer a truncated model on [x0,1]
use_truncated = true;       % set false to use untruncated shifted exponential

%% ---------- 4) Prepare x-grid and model curve ----------
x_grid = linspace(min(x), max(x), 400);

if use_truncated
    % Truncated shifted exponential on [x0, 1]
    xmax = 1;
    Z = 1 - exp(-lambda*(xmax - x0)); % normalization over [x0, xmax]
    model_pdf = zeros(size(x_grid));
    mask = (x_grid >= x0) & (x_grid <= xmax);
    model_pdf(mask) = (lambda * exp(-lambda*(x_grid(mask)-x0))) / Z;
else
    % Untruncated shifted exponential on [x0, inf)
    model_pdf = zeros(size(x_grid));
    mask = (x_grid >= x0);
    model_pdf(mask) = lambda * exp(-lambda*(x_grid(mask)-x0));
end

%% ---------- 5) Plot histogram + overlay ----------
figure('Position',[100 100 720 520]); hold on;

if use_pdf
    h = histogram(x, nbins, 'Normalization','pdf', 'DisplayStyle','bar');
    y_label = 'PDF';
    overlay_y = model_pdf;
else
    h = histogram(x, nbins, 'Normalization','count', 'DisplayStyle','bar');
    y_label = 'Frequency';
    % Scale PDF to counts: counts ≈ pdf * N * binWidth
    binWidth = h.BinWidth;
    overlay_y = model_pdf * numel(x) * binWidth;
end

% Overlay
plot(x_grid, overlay_y, 'LineWidth', 2);

% Cosmetics
xlabel('Granule fraction');
ylabel(y_label);
title('Granule fraction distribution with minimal-model overlay');

if log_y
    set(gca, 'YScale', 'log');
end

% Show fitted params
if use_truncated
    model_name = 'shifted exp (truncated to [0,1])';
else
    model_name = 'shifted exp (untruncated)';
end

legend({'Data', sprintf('%s: \\lambda=%.3g', model_name, lambda)}, ...
       'Location','best');

grid on;

%% ---------- 6) Optional: also show CCDF on a separate figure (often cleaner on log scale) ----------
% CCDF (survival) is often visually more stable than histogram on log-y
do_ccdf = true;
if do_ccdf
    figure('Position',[860 100 720 520]); hold on;
    % Empirical survival
    x_sorted = sort(x);
    N = numel(x_sorted);
    surv = (N:-1:1)'/N; % P(X >= x_sorted(i))
    semilogy(x_sorted, surv, '.', 'MarkerSize', 10);

    % Model survival
    if use_truncated
        xmax = 1;
        Z = 1 - exp(-lambda*(xmax - x0));
        S_model = ones(size(x_grid));
        % For x>=x0: S(x)= (exp(-lambda(x-x0)) - exp(-lambda(xmax-x0))) / Z
        S_model = zeros(size(x_grid));
        mask = (x_grid < x0);
        S_model(mask) = 1;
        mask = (x_grid >= x0) & (x_grid <= xmax);
        S_model(mask) = (exp(-lambda*(x_grid(mask)-x0)) - exp(-lambda*(xmax-x0))) / Z;
        mask = (x_grid > xmax);
        S_model(mask) = 0;
    else
        S_model = ones(size(x_grid));
        mask = (x_grid >= x0);
        S_model(mask) = exp(-lambda*(x_grid(mask)-x0));
        S_model(x_grid < x0) = 1;
    end
    semilogy(x_grid, S_model, 'LineWidth', 2);

    xlabel('Granule fraction');
    ylabel('P(Fraction \ge x)');
    title('CCDF (survival) with minimal-model overlay');
    grid on;
    legend({'Empirical CCDF', sprintf('Model CCDF (x0=%.3g, \\lambda=%.3g)', x0, lambda)}, ...
        'Location','best');
end
