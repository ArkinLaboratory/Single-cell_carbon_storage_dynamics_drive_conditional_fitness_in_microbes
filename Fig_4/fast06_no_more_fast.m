row_data = total_cell_data_collect(14, :); 

figure;
histogram(row_data, 'BinWidth', 10);

title('histogram gpr');
xlabel('Granule Prod Rate');
ylabel('Freq');

grid on;
%% 
indices = find(total_cell_data_collect(14, :) > 33);

fast_index = unique(total_cell_data_collect(10, indices));
fast_count = length(fast_index);

fprintf('fast count (threshold 33): %f\n', fast_count);

%% 

normal_time_collect = [];
fast_normal_fast_length = length(fast_normal_fast_collect(1,:));
for i = 1:length(fast_normal_fast)
% for i = 1
    normal_time = 0;
    cell_index = fast_normal_fast(i);
    cell_cycles = find(fast_normal_fast_collect(10, :) == cell_index);
    trajec_length = length(cell_cycles);
    for j = 1:trajec_length
        cycle_num = cell_cycles(j);

        next_cycle = cycle_num + 1;
        if next_cycle < fast_normal_fast_length && fast_normal_fast_collect(10, next_cycle) == fast_normal_fast(i) && fast_normal_fast_collect(7, cycle_num) > 20
            % fprintf('%f\n', fast_normal_fast(i));
            while cycle_num < fast_normal_fast_length && fast_normal_fast_collect(7, cycle_num + 1)< 20 && fast_normal_fast_collect(10, cycle_num + 1) == fast_normal_fast_collect(10, cycle_num)
                % fprintf('%f\n', fast_normal_fast_collect(6, cycle_num + 1));
                normal_time = normal_time + fast_normal_fast_collect(6, cycle_num + 1);
                cycle_num = cycle_num + 1;
            end
        end
    end      
    % fprintf('%f\n',normal_time);
    normal_time_collect = [normal_time_collect; cell_index, normal_time];
end

%% 
figure;
hold on

bin_length = 5;
row = 12;

values = total_cell_data_collect(row, :);
% values = total_cell_data_collect(row, :);

h2 = histogram(values, 'BinWidth', bin_length, 'Normalization', 'probability');
% h201 = histogram(values_end, 'BinWidth', bin_length, 'Normalization', 'probability');
% h2 = histogram(values, 'BinWidth', bin_length);


% legend([h2, h201], 'intial granule size', 'end granule size');

% switch time
ylims = ylim;
line([276 276], ylims, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
line([496 496], ylims, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 2);
line([710 710], ylims, 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2);
legend('# of fast', 't = 276', 't = 496', 't = 710');
% % legend('# of all cells', 't = 253', 't = 463', 't = 692');

% line([259 259], ylims, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
% legend('# of all cells', 't = 259')
% legend([h1, h2], 'Total Cell', 'Fast');

% title('fast N conc 0.1 to 0.001 (Threshold 33)');
title('C starv Lane 2');
% title('4-20 all cells');

xlabel('time');
% xlabel('% Concentration of N');

% ylabel('Percentage');
ylabel('Num');

grid on;
hold off

%% 
lane1_values = values;
%% 
lane2_values = values;
%% 
figure;
h2 = histogram(lane1_values, 'BinWidth', bin_length, 'Normalization', 'probability');
hold on;
h201 = histogram(lane2_values, 'BinWidth', bin_length, 'Normalization', 'probability');



% h2 = histogram(values, 'BinWidth', bin_length);


% switch time
ylims = ylim;
line([276 276], ylims, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
line([496 496], ylims, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 2);
line([710 710], ylims, 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2);
legend('lane1', 'lane2', 't = 276', 't = 496', 't = 710');
% % legend('# of all cells', 't = 253', 't = 463', 't = 692');

% line([259 259], ylims, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
% legend('# of all cells', 't = 259')
% legend([h1, h2], 'Total Cell', 'Fast');

% title('fast N conc 0.1 to 0.001 (Threshold 33)');
title('C starv both');
% title('4-20 all cells');

xlabel('time');
% xlabel('% Concentration of N');

% ylabel('Percentage');
ylabel('Num');

grid on;
hold off

%% 
bucket01 = [];
i = 0;
step = 10;
while i < 851
    bucket01 = [bucket01, i];
    i = i + step;
end

ini_count = sum(total_cell_data_collect(11, :) < step);

bucket01(2,1) = ini_count;

for j = 2:length(bucket01)
    i = bucket01(1,j);
    start_count = sum(total_cell_data_collect(11, :) < i & total_cell_data_collect(11, :) > (i-step));
    end_count = sum(total_cell_data_collect(12, :) < i & total_cell_data_collect(12, :) > (i-step));
    total_count = bucket01(2,j-1) + start_count - end_count;

    bucket01(2,j) = total_count;
    bucket01(3,j) = start_count;
    bucket01(4,j) = end_count;
    bucket01(5,j) = start_count / total_count;
    bucket01(6,j) = end_count / total_count;
end
%% 
lane1_bucket = bucket01;
%% 
lane2_bucket = bucket01;
%% 
figure;
end_time = length(bucket01)-3;
plot(bucket01(1,1:end_time), bucket01(5,1:end_time))
hold on;
plot(bucket01(1,1:end_time), bucket01(6,1:end_time))
legend('start count', 'end count')
title('lane2')
xlabel('time')
ylabel('count over total cells')
%% 
figure;
end_time = length(bucket01)-3;
row = 6;
bar(lane1_bucket(1, 1:end_time), lane1_bucket(row, 1:end_time)./step, 'BarWidth', 1, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on;
bar(lane2_bucket(1, 1:end_time), lane2_bucket(row, 1:end_time)./step, 'BarWidth', 1, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

ylims = ylim;

line([0 0], ylims, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 2);
line([243 243], ylims, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
line([464 464], ylims, 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2);

% line([0 0], ylims, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 2);
% line([272 272], ylims, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
% line([488 488], ylims, 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2);

legend('WT', 'ΔPhaC', '2C2N', '0C2N', '2e-4C 2e-4N');
title('C starv end')
xlabel('time')
ylabel('count over total cells')
%% 
N_aver_lane2 = mean(total_cell_data_collect(6,:));
%% 
figure;
histogram(total_cell_data_collect(6,:),'BinWidth',1)
%% 
N_lane1_bucket = bucket01;
%% 
N_lane2_bucket = bucket01;
%% 
figure;
end_time = length(bucket01)-3;
row = 6;
bar(N_lane1_bucket(1, 1:end_time), N_lane1_bucket(row, 1:end_time), 'BarWidth', 1);
hold on;
bar(N_lane2_bucket(1, 1:end_time), N_lane2_bucket(row, 1:end_time), 'BarWidth', 1);
ylims = ylim;
line([0 0], ylims, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
line([272 272], ylims, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 2);
line([488 488], ylims, 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2);
legend('lane1', 'lane2', '2C2N', '2C0N', '2e-3 CN');
title('C starv end')
xlabel('time')
ylabel('count over total cells')
%% 
row10 = total_cell_data_collect(10, :);

[unique_values, ~, idx] = unique(row10, 'stable');
last_occurrences = arrayfun(@(x) find(row10 == x, 1, 'last'), unique_values);

total_cell_data_collect(:, last_occurrences) = []; 