%% Minimal granule segregation + growth model (faster, subplot)

rng(1);

%% Parameters
params.tau_mean   = 90;       % mean cell-cycle time τ [min]
params.tau_cv     = 0.10;     % CV of τ
params.p_keep     = 0.5;      % prob granule stays in tracked lineage each division
params.s0         = 1.0;      % size at acquisition
params.g          = 0.03;     % linear growth rate [size/min]
params.growthLaw  = "linear";    % <<< DEFAULT NOW "exp"  ("linear" or "exp")
params.gamma      = 0.01;     % exp growth rate [/min] if growthLaw="exp"

% Acquisition delay (lognormal): median ~ 40 min
params.Tacq_median = 40;
params.Tacq_sigma  = 0.35;
params.Tacq_mu     = log(params.Tacq_median);

%% Simulation scale (FAST defaults)
sim.N_lineages   = 800;
sim.T_total      = 8000;
sim.T_burnin     = 1500;
sim.sampleEvery  = 50;

out = simulate_lineages_fast(params, sim);

%% ---- Build one figure with subplots ----
figure('Position',[100 100 1200 700]);

% (5) Size PDF at steady state
subplot(1,2,1);
edgesS = linspace(min(out.S_sample), prctile(out.S_sample, 99.5), 60);
h = histogram(out.S_sample, edgesS, 'Normalization', 'pdf');
xlabel('s (a.u.)'); ylabel('PDF');
title('Granule fraction distribution');

% <<< Make PDF y-axis log-scaled >>>
set(gca, 'YScale', 'log');
grid on;

% Avoid log(0): set a small positive lower limit
% pick something below smallest positive bin height, but not too tiny
binHeights = h.Values;
minPos = min(binHeights(binHeights > 0));
if isempty(minPos)
    minPos = 1e-6;
end
ylim([minPos/2, max(binHeights)*1.2]);

% (6) Survival plot + fit depending on growth law
subplot(1,2,2);

s0 = params.s0;
sFit = out.S_sample(out.S_sample >= s0);

% Empirical survival function on a grid
sGrid = linspace(s0, prctile(sFit, 99.5), 80);
survS = arrayfun(@(s) mean(sFit >= s), sGrid);

if params.growthLaw == "linear"
    % ---- shifted exponential tail in (s - s0) ----
    x = (sFit - s0);
    x = x(x > 0);
    lambda_hat = 1 / mean(x);  % MLE for Exp on x

    semilogy(sGrid, survS, 'o-'); hold on; grid on;
    semilogy(sGrid, exp(-lambda_hat*(sGrid - s0)), '--', 'LineWidth', 2);
    xlabel('s'); ylabel('P(S \ge s)');
    title('Survival: shifted exponential (linear growth)');
    legend('Empirical', sprintf('fit exp(-\\lambda(s-s0)), \\lambda=%.3f', lambda_hat), ...
        'Location','southwest');

elseif params.growthLaw == "exp"
    % ---- power-law (Pareto) tail in s ----
    % Model: P(S >= s) = (s0/s)^alpha for s >= s0
    % MLE for alpha (Pareto Type I with xmin=s0):
    % alpha_hat = n / sum(log(s_i / s0))
    ratio = sFit / s0;
    alpha_hat = numel(sFit) / sum(log(ratio));

    loglog(sGrid, survS, 'o-'); hold on; grid on;
    loglog(sGrid, (s0 ./ sGrid) .^ alpha_hat, '--', 'LineWidth', 2);

    xlabel('s'); ylabel('P(S \ge s)');
    title('Survival: power-law tail (exp growth)');
    legend('Empirical', sprintf('fit (s0/s)^\\alpha, \\alpha=%.3f', alpha_hat), ...
        'Location','southwest');

else
    text(0.1,0.5,'Unknown growthLaw','Units','normalized');
    axis off;
end

sgtitle('Minimal segregation + growth model (single-lineage tracking)');

%% Export to PDF (optional)
% Change filename/path as you like:
exportgraphics(gcf, 'S5_minimal_model.pdf', 'ContentType','vector');

%% Print summary
fprintf('--- Summary ---\n');
fprintf('N_lineages=%d, T_total=%.0f min, sampleEvery=%.0f min\n', sim.N_lineages, sim.T_total, sim.sampleEvery);
fprintf('p_keep=%.2f, tau_mean=%.1f min, Tacq_median=%.1f min\n', params.p_keep, params.tau_mean, params.Tacq_median);
fprintf('Median lifetime: K=%.1f gen, T=%.1f min\n', median(out.K_life), median(out.T_life));
fprintf('Steady-state sampled granule-bearing fraction ≈ %.3f\n', mean(out.hasGranule_sample));

%% ---------------- Local functions ----------------
function out = simulate_lineages_fast(params, sim)
    N = sim.N_lineages;

    hasGranule = false(N,1);
    age = zeros(N,1);
    K = zeros(N,1);
    timeToAcquire = draw_Tacq(params, N);
    nextDiv = draw_tau(params, N);

    K_life = zeros(0,1);
    T_life = zeros(0,1);
    S_sample = zeros(0,1);
    hasGranule_sample = false(0,1);

    t = 0;
    nextSample = sim.T_burnin;

    while t < sim.T_total
        tNextDiv = min(nextDiv);
        tNext = min(tNextDiv, nextSample);

        dt = tNext - t;

        if dt > 0
            age(hasGranule) = age(hasGranule) + dt;
            timeToAcquire(~hasGranule) = timeToAcquire(~hasGranule) - dt;
        end

        newly = (~hasGranule) & (timeToAcquire <= 0);
        if any(newly)
            hasGranule(newly) = true;
            age(newly) = 0;
            K(newly) = 0;
            timeToAcquire(newly) = NaN;
        end

        if abs(tNext - nextSample) < 1e-9
            s = size_from_age(params, age(hasGranule));
            S_sample = [S_sample; s]; %#ok<AGROW>
            hasGranule_sample = [hasGranule_sample; hasGranule]; %#ok<AGROW>
            nextSample = nextSample + sim.sampleEvery;
        end

        divIdx = find(abs(nextDiv - tNext) < 1e-9);
        if ~isempty(divIdx)
            idx = divIdx(:);

            gIdx = idx(hasGranule(idx));
            if ~isempty(gIdx)
                K(gIdx) = K(gIdx) + 1;

                keepMask = rand(size(gIdx)) <= params.p_keep;
                loseIdx = gIdx(~keepMask);

                if ~isempty(loseIdx)
                    K_life = [K_life; K(loseIdx)]; %#ok<AGROW>
                    T_life = [T_life; age(loseIdx)]; %#ok<AGROW>

                    hasGranule(loseIdx) = false;
                    age(loseIdx) = 0;
                    K(loseIdx) = 0;
                    timeToAcquire(loseIdx) = draw_Tacq(params, numel(loseIdx));
                end
            end

            nextDiv(idx) = tNext + draw_tau(params, numel(idx));
        end

        t = tNext;
    end

    out.K_life = K_life;
    out.T_life = T_life;
    out.S_sample = S_sample;
    out.hasGranule_sample = hasGranule_sample;
end

function tau = draw_tau(params, n)
    if params.tau_cv <= 0
        tau = params.tau_mean * ones(n,1);
        return;
    end
    cv = params.tau_cv;
    mu = log(params.tau_mean / sqrt(1 + cv^2));
    sigma = sqrt(log(1 + cv^2));
    tau = lognrnd(mu, sigma, [n,1]);
end

function Tacq = draw_Tacq(params, n)
    Tacq = lognrnd(params.Tacq_mu, params.Tacq_sigma, [n,1]);
end

function s = size_from_age(params, age)
    if params.growthLaw == "linear"
        s = params.s0 + params.g * age;
    elseif params.growthLaw == "exp"
        s = params.s0 * exp(params.gamma * age);
    else
        error('growthLaw must be "linear" or "exp"');
    end
end
