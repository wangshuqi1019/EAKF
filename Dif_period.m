clear;
addpath('data_shenzhen/')

[num1,~,~] = xlsread('2017-2019.xlsx');
%[num2,~,~] = xlsread('2023.xlsx');
stateNums = zeros(56*2,1); % 一共2017-2019年跨越两年的数据，每年56周（为了避免误差多考虑4周）
temp1 = num1(:,5);
% temp2 = num2(:,5);
stateNums(1:length(temp1)) = temp1;

% % ILI+
% Observation = stateNums(:,1)'*0.1;
% Observation = Observation';

%% test
iter = 0;
mu = 0.1;
init = 0.65;
param_variance = 0.005;
betavar = 0.001;
betamean = 0.02;
oev = 0.2;

% load('result/S0_tri.mat');
rng(300);
% 生成 LHS 样本
n_samples = 100;
lhs_samples = lhsdesign(n_samples, 2);
% 将 LHS 样本映射到 [0.65, 0.75] 范围内
s1_samples = 0.65 + lhs_samples(:, 1) * (0.75 - 0.65);
s2_samples = 0.65 + lhs_samples(:, 2) * (0.75 - 0.65);
S0_tri = zeros(100,2);
for i = 1:n_samples
    s1 = s1_samples(i);
    s2 = s2_samples(i);
    S0_tri(i,1) = s1;
    S0_tri(i,2) = s2;   
end
save('result/S0_tri.mat','S0_tri');

% flag = Flu_dif_period_adj(stateNums, mu, S0_tri, param_variance, betamean, betavar);
%% 1
[MSE_list1, result_save1, result_w1, result_beta1] = ...
  Flu_dif_period1(stateNums, mu, S0_tri(1:50,:), param_variance, betamean, betavar, oev);
mse = sum(MSE_list1);
disp(['var = ' num2str(param_variance) ', betamean = ' num2str(betamean) ', betavar = ' num2str(betavar) ', oev = ' num2str(oev) ', mse = ' num2str(mse)]);
save('result/result_w1.mat','result_w1');
save('result/result_beta1.mat','result_beta1');

mean_result_save1 = {};
for i = 1:length(result_save1)
    % 获取 1x2 cell 数组中的两个 29x10000x56 数组
    data1 = result_save1{i}{1};
    data2 = result_save1{i}{2};
    
    % 计算均值和标准差
    mean_data1 = mean(data1, 2);  % 对第二个维度求均值
    mean_data1 = squeeze(mean_data1);
    mean_data2 = mean(data2, 2);  % 对第二个维度求均值
    mean_data2 = squeeze(mean_data2);
    
    % 存储计算结果
    mean_result_save1{i} = {mean_data1, mean_data2};
end
save('result/mean_result_save1.mat','mean_result_save1')
clear clear MSE_list1
clear -regexp ^result
clear -regexp ^mean
clear -regexp ^data

%% 2
[MSE_list2, result_save2, result_w2, result_beta2] = ...
  Flu_dif_period1(stateNums, mu, S0_tri(51:100,:), param_variance, betamean, betavar, oev);
mse = sum(MSE_list2);
disp(['var = ' num2str(param_variance) ', betamean = ' num2str(betamean) ', betavar = ' num2str(betavar) ', oev = ' num2str(oev) ', mse = ' num2str(mse)]);
save('result/result_w2.mat','result_w2');
save('result/result_beta2.mat','result_beta2');

mean_result_save2 = {};
for i = 1:length(result_save2)
    % 获取 1x2 cell 数组中的两个 29x10000x56 数组
    data1 = result_save2{i}{1};
    data2 = result_save2{i}{2};
    
    % 计算均值和标准差
    mean_data1 = mean(data1, 2);  % 对第二个维度求均值
    mean_data1 = squeeze(mean_data1);
    mean_data2 = mean(data2, 2);  % 对第二个维度求均值
    mean_data2 = squeeze(mean_data2);
    
    % 存储计算结果
    mean_result_save2{i} = {mean_data1, mean_data2};
end
save('result/mean_result_save2.mat','mean_result_save2')

%% 保存结果
load('result/mean_result_save1.mat');
load('result/mean_result_save2.mat');
mean_forecast = {};
for i=1:50
    mean_forecast{i} = mean_result_save1{i};
end
for i=51:100
    mean_forecast{i} = mean_result_save2{i-50};
end

%% 计算均值和方差
% 初始化结果存储变量
state_index = 29;
num_repeats = 100;
num_periods = 2;
num_weeks = 56;

% 存储两个时间周期的均值和方差
mean_values = zeros(num_periods, num_weeks);
std_values = zeros(num_periods, num_weeks);

% 遍历每个时间周期
for period = 1:num_periods
    % 提取当前时间周期的所有实验数据
    all_experiments = zeros(num_repeats, num_weeks);
    for repeat = 1:num_repeats
        all_experiments(repeat, :) = mean_forecast{repeat}{period}(state_index, :);
    end
    
    % 计算均值和标准差
    mean_values(period, :) = mean(all_experiments, 1);
    std_values(period, :) = std(all_experiments, 0, 1);
    
end
mean_values = mean_values';
std_values = std_values';
save('result/mean_values.mat','mean_values');
save('result/std_values.mat','std_values');

%% 计算均值和置信区间
% 初始化结果存储变量
state_index = 29;
num_repeats = 100;
num_periods = 2;
num_weeks = 56;

% 存储两个时间周期的均值和置信区间
mean_values = zeros(num_weeks, num_periods);
CI_lower_values = zeros(num_weeks, num_periods);
CI_upper_values = zeros(num_weeks, num_periods);

% 遍历每个时间周期
for period = 1:num_periods
    % 提取当前时间周期的所有实验数据
    all_experiments = zeros(num_repeats, num_weeks);
    for repeat = 1:num_repeats
        all_experiments(repeat, :) = mean_forecast{repeat}{period}(state_index, :);
    end
    
    % 计算均值和置信区间
    for week = 1:num_weeks
        data = all_experiments(:, week);
        SEM = std(data) / sqrt(length(data)); % 计算标准误差
        ts = tinv([0.025, 0.975], length(data) - 1); % 计算T-Score
        CI = mean(data) + ts * SEM; % 计算置信区间
        mean_values(week, period) = mean(data);
        std_values(week, period) = std(data);
        CI_lower_values(week, period) = CI(1);
        CI_upper_values(week, period) = CI(2);
    end
end
save('result/mean_values.mat','mean_values');
save('result/CI_lower_values.mat','CI_lower_values');
save('result/CI_upper_values.mat','CI_upper_values');

