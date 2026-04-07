%% CN vs PA pulsed-regime model + analytic overlays (C-limit + N-limit)
clear; clc;

% =======================
% Parameters
% =======================
r_CN   = 0.2773;          % CN growth rate in abundance
T_end  = 100;             % total simulation time
dt     = 0.1;             % time step
t_vec  = 0:dt:T_end;

duties = 0:0.1:1;         % duty cycle D grid
x_vals = 0:10:100;        % PA speed advantage (%)
n_x    = numel(x_vals);

T_gen_CN = log(2)/r_CN;
fprintf('CN division time = %.3f\n', T_gen_CN);

log_C_all   = cell(n_x,1);
log_N_all   = cell(n_x,1);
periods_all = cell(n_x,1);

% =======================
% Simulation
% =======================
for xi = 1:n_x
    x    = x_vals(xi);
    r_PA = r_CN * (1 + x/100);
    T_gen_PA = log(2)/r_PA;

    minP = ceil( max(T_gen_CN, T_gen_PA) / 10 ) * 10;
    fprintf('  x = %3d%%, PA division time = %.3f, minP = %2d\n', x, T_gen_PA, minP);

    periods_x = minP:10:100;
    nP        = numel(periods_x);
    periods_all{xi} = periods_x;

    log_C = nan(nP, numel(duties));
    log_N = nan(nP, numel(duties));

    % half-rate phase duration in your code (named t_two_*_half)
    t_two_CN_half = 4*log(2)/r_CN;
    t_two_PA_half = 4*log(2)/r_PA;

    % -------- C-limit simulation --------
    for ip = 1:nP
        p = periods_x(ip);
        for id = 1:numel(duties)
            d  = duties(id);
            R  = double(mod(t_vec, p) < d*p);   % 1 in abundance, 0 in starvation
            dR = [0, diff(R)];

            bCN = ones(size(t_vec));
            bPA = ones(size(t_vec));
            ph_CN = 0;
            ph_PA = 0;

            for i = 2:numel(t_vec)
                % end of abundance: CN gets a fold change 1.3 (one event per cycle)
                if dR(i) == -1
                    bCN(i-1) = bCN(i-1) * 1.3;
                end
                % start of abundance: reset "ph" timers
                if dR(i) == +1
                    ph_CN = t_two_CN_half;
                    ph_PA = t_two_PA_half;
                end

                % CN update
                if ph_CN > 0
                    bCN(i) = bCN(i-1) * exp((r_CN/2)*dt);
                    ph_CN = ph_CN - dt;
                elseif R(i)
                    bCN(i) = bCN(i-1) * exp(r_CN*dt);
                else
                    bCN(i) = bCN(i-1);
                end

                % PA update
                if ph_PA > 0
                    bPA(i) = bPA(i-1) * exp((r_PA/2)*dt);
                    ph_PA = ph_PA - dt;
                elseif R(i)
                    bPA(i) = bPA(i-1) * exp(r_PA*dt);
                else
                    bPA(i) = bPA(i-1);
                end
            end

            log_C(ip,id) = log(bCN(end)/bPA(end));  % natural log
        end
    end

    % -------- N-limit simulation --------
    for ip = 1:nP
        p = periods_x(ip);
        for id = 1:numel(duties)
            d  = duties(id);
            R  = double(mod(t_vec, p) < d*p);
            dR = [0, diff(R)];

            bCN = ones(size(t_vec));
            bPA = ones(size(t_vec));

            wait_PA      = false;
            elapsed_PA   = 0;
            gen_timer_PA = 0;

            for i = 2:numel(t_vec)
                % start of abundance: PA enters waiting state for one generation
                if dR(i) == +1
                    wait_PA      = true;
                    gen_timer_PA = log(2)/r_PA;
                    elapsed_PA   = 0;
                end

                % PA update
                if R(i)
                    if wait_PA
                        bPA(i)     = bPA(i-1);
                        elapsed_PA = elapsed_PA + dt;
                        if elapsed_PA >= gen_timer_PA
                            wait_PA = false;
                        end
                    else
                        bPA(i) = bPA(i-1) * exp(r_PA*dt);
                    end
                else
                    bPA(i) = bPA(i-1);
                end

                % CN update
                if R(i)
                    bCN(i) = bCN(i-1) * exp(r_CN*dt);
                else
                    bCN(i) = bCN(i-1);
                end
            end

            log_N(ip,id) = log(bCN(end)/bPA(end));  % natural log
        end
    end

    log_C_all{xi} = log_C;
    log_N_all{xi} = log_N;
end

% =======================
% Plotting
% =======================
x_plot = [1,3,5,7,9];  % choose 5 x indices
figure('Position',[100 100 200*numel(x_plot) 500]);
tl = tiledlayout(2,numel(x_plot),'TileSpacing','compact','Padding','compact');

% color axis centered at 0 across both rows
allC = cellfun(@(M) M(:), log_C_all(x_plot),'UniformOutput',false);
allN = cellfun(@(M) M(:), log_N_all(x_plot),'UniformOutput',false);
v    = vertcat(allC{:}, allN{:});
L    = max(abs(v));
cax  = [-L, L];

colormap(greenZeroDiverging(256));

for k = 1:numel(x_plot)
    xi = x_plot(k);
    P  = periods_all{xi};
    Cm = log_C_all{xi};
    Nm = log_N_all{xi};

    % ---------- Common analytic helpers for this xi ----------
    rPA = r_CN * (1 + x_vals(xi)/100);

    duty_fine = linspace(0.01, 0.99, 500); % avoid D=0 and D=1
    Pmin = min(P);
    Pmax = max(P);

    % ===== N-limit analytic: kink + 0-line =====
    L_lag = log(2)/rPA;                 % kink: DP = L_lag
    P_kink_N = L_lag ./ duty_fine;

    den = (rPA - r_CN);
    if den > 1e-12
        DP0_N = log(2) / den;           % 0-line: DP = ln2/(rPA - rCN)
        P_zero_N = DP0_N ./ duty_fine;
    else
        P_zero_N = nan(size(duty_fine)); % x=0 => no nontrivial 0-line
    end

    P_kink_N(P_kink_N < Pmin | P_kink_N > Pmax) = NaN;
    P_zero_N(P_zero_N < Pmin | P_zero_N > Pmax) = NaN;

    % ===== C-limit analytic: 0-line + two "ph" boundaries =====
    % Your ph duration:
    H_CN = 4*log(2)/r_CN;
    H_PA = 4*log(2)/rPA;

    % C-limit DP0 has piecewise candidate forms (matching your ph_* + 1.3 mechanism)
    DP0_C = NaN;

    % Candidate A: H_PA < DP0 <= H_CN
    DP0_A = log(1.3*16) / rPA;              % 24 = 1.5 * 16 (since exp((rPA/2)*H_PA)=16)
    if (DP0_A > H_PA) && (DP0_A <= H_CN)
        DP0_C = DP0_A;
    else
        % Candidate B: DP0 > H_CN
        if den > 1e-12
            DP0_B = log(1.5) / den;
            if DP0_B > H_CN
                DP0_C = DP0_B;
            end
        end
    end

    P_zero_C = (DP0_C ./ duty_fine);

    % "ph" boundaries (mechanism switch lines)
    P_kink_PA_C = H_PA ./ duty_fine;    % DP = H_PA
    P_kink_CN_C = H_CN ./ duty_fine;    % DP = H_CN

    P_zero_C(P_zero_C < Pmin | P_zero_C > Pmax) = NaN;
    P_kink_PA_C(P_kink_PA_C < Pmin | P_kink_PA_C > Pmax) = NaN;
    P_kink_CN_C(P_kink_CN_C < Pmin | P_kink_CN_C > Pmax) = NaN;

    % ---------- C-limit tile ----------
    nexttile(k);
    imagesc(duties, P, Cm);
    axis xy; caxis(cax);
    xlabel('Duty'); if k==1, ylabel('Period'); end
    title(sprintf('C-limit, x = %d%%', x_vals(xi)));
    hold on;
    % contour(duties, P, Cm, [0 0], 'k-', 'LineWidth', 0.9); % numeric 0-line

    % C-limit analytic overlays (pink)
    plot(duty_fine, P_zero_C,    'r-',  'LineWidth', 1.4); % analytic 0-line for C-limit
    % plot(duty_fine, P_kink_PA_C, 'm--', 'LineWidth', 1.0); % DP=H_PA (optional)
    % plot(duty_fine, P_kink_CN_C, 'm--', 'LineWidth', 1.0); % DP=H_CN (optional)
    hold off;

    % ---------- N-limit tile ----------
    nexttile(numel(x_plot)+k);
    imagesc(duties, P, Nm);
    axis xy; caxis(cax);
    xlabel('Duty'); if k==1, ylabel('Period'); end
    title(sprintf('N-limit, x = %d%%', x_vals(xi)));
    hold on;
    % contour(duties, P, Nm, [0 0], 'k-', 'LineWidth', 0.9); % numeric 0-line

    % N-limit analytic overlays (white)
    plot(duty_fine, P_zero_N, 'r-',  'LineWidth', 1.2);    % analytic 0-line for N-limit
    % plot(duty_fine, P_kink_N, 'w--', 'LineWidth', 1.0);    % kink DP=L (optional)
    hold off;
end

% Shared colorbar
try
    cb = colorbar(tl,'eastoutside');
    cb.Label.String = 'log(CN/PA)  (natural log)';
catch
    nexttile(numel(x_plot));
    colorbar;
end

sgtitle('Comparison of average growth ratios at different speed advantages');

% =======================
% Colormap helper (green at zero)
% =======================
function cmap = greenZeroDiverging(n)
    if nargin<1, n = 256; end
    half = floor(n/2);
    % negative: blue -> green
    neg = [ zeros(half,1), linspace(0,1,half)', linspace(1,0,half)' ];
    % positive: green -> yellow
    pos = [ linspace(0,1,n-half)', ones(n-half,1), zeros(n-half,1) ];
    cmap = [neg; pos];
end
