figure;
h = histogram(total_cell_data_collect(14, :), 'BinWidth', 0.1);
% h = h(h ~= 0);
% h = histogram(total_cell_data_collect(16, :), 'BinWidth', 0.005);
set(gca, 'YScale', 'log');

h.EdgeColor = h.FaceColor;
xlim([-50 100])
xlabel('Granule Pridcution Rate');
ylabel('Frequency (log scale)');
% title('Histogram with Logarithmic Scale on Y-axis');