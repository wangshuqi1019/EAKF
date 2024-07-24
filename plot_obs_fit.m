function plot_obs_fit(Observation, mean_values, std_values)
    % 分割数据为每年的数据
    obs_year1 = Observation(1:56);
    obs_year2 = Observation(57:end);

    mean_year1 = mean_values(:, 1);
    mean_year2 = mean_values(:, 2);

    std_year1 = std_values(:, 1);
    std_year2 = std_values(:, 2);

    % 清除画布
    clf;

    % 创建图形和子图
    figure;

    % 绘制第一年的数据
    subplot(2, 1, 1);
    hold on;
    plot(1:56, obs_year1, 'o-', 'DisplayName', 'Observation Year 1');
    plot(1:56, mean_year1, 'r-', 'DisplayName', 'Fit Mean Year 1');
    fill([1:56, fliplr(1:56)], [mean_year1 - std_year1; flipud(mean_year1 + std_year1)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Fit Std Year 1');
    ylabel('Observation and Fit Data');
    legend;
    title('Year 1: Observation vs. Fit');
    hold off;

    % 绘制第二年的数据
    subplot(2, 1, 2);
    hold on;
    plot(1:56, obs_year2, 'o-', 'DisplayName', 'Observation Year 2');
    plot(1:56, mean_year2, 'r-', 'DisplayName', 'Fit Mean Year 2');
    fill([1:56, fliplr(1:56)], [mean_year2 - std_year2; flipud(mean_year2 + std_year2)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Fit Std Year 2');
    ylabel('Observation and Fit Data');
    xlabel('Weeks');
    legend;
    title('Year 2: Observation vs. Fit');
    hold off;
end
