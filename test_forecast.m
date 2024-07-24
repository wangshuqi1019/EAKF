clear;
addpath('data_shenzhen/')

%% fitting
%% load fitting data
nperiod = 2;
nweek = 56;
[num1,~,~] = xlsread('2017-2019.xlsx');
stateNums = zeros(nweek*nperiod,1); % 56 weeks, 2 years
temp1 = num1(:,5);
stateNums(1:length(temp1)) = temp1;

%% fitting test
nSim = 50;  % origin 100, repeating times, limited by server
ns = 1000;  % number of particles

mu = 0.1;  % ILI+ mapping
init = 0.65;  % init S0
param_variance = 0.01;
betavar = 0.001;
betamean = 0.01;
oev = 1e-5;
% 直接预测beta，没有用R0计算beta
load('result/S0_tri.mat');
[MSE_list, result_save, result_w, result_beta] = ...
    Flu_dif_period_test(stateNums, mu, S0_tri(1:nSim,:), param_variance, betamean, betavar, oev, ns);
mse = sum(MSE_list);
disp(mse);
% disp(['var = ' num2str(param_variance) ', betamean = ' num2str(betamean) ', betavar = ' num2str(betavar) ', oev = ' num2str(oev) ', mse = ' num2str(mse)]);

%% save fitting results
para.param_variance=param_variance;
para.mu=mu;
para.betamean=betamean;
para.betavar=betavar;
para.oev=oev;
para.init=init;
save('forecast/para.mat','para');

mean_result_save = {};
for i = 1:length(result_save)  % size of result_save: 1*50 cell
    data1 = result_save{i}{1};  % size of data1: 32x10000x56
    data2 = result_save{i}{2};  % size of data2: 32x10000x56
    mean_data1 = mean(data1, 2);  % 对第二个维度求均值
    mean_data1 = squeeze(mean_data1);
    mean_data2 = mean(data2, 2);  % 对第二个维度求均值
    mean_data2 = squeeze(mean_data2);
    mean_result_save{i} = {mean_data1, mean_data2};
end

mean_forecast = {};
for i=1:nSim
    mean_forecast{i} = mean_result_save{i};
end
mean_values = zeros(nweek, nperiod);
std_values = zeros(nweek, nperiod);
for period = 1:nperiod
    all_experiments = zeros(nSim, nweek);
    for repeat = 1:nSim
        all_experiments(repeat, :) = mean_forecast{repeat}{period}(29, :);  % symptomatic incidence
    end
    mean_values(:, period) = mean(all_experiments, 1);
    std_values(:, period) = std(all_experiments, 0, 1);
end
save('forecast/mean_values.mat','mean_values');
save('forecast/std_values.mat','std_values');
% save('forecast/result_w.mat','result_w');
% save('forecast/result_beta.mat','result_beta');

%% forecasting
%% load forecasting data
w1 = zeros(nSim, ns, nweek);
for iter = 1:nSim
    season1_data = result_w{iter}{1};
    for week = 1:nweek
        w1(iter, :, week) = season1_data(:, week);
    end
end
% save('result/w1.mat','w1');
w2 = zeros(nSim, ns, nweek);
for iter = 1:nSim
    season2_data = result_w{iter}{2};
    for week = 1:nweek
        w2(iter, :, week) = season2_data(:, week);
    end
end
% save('result/w2.mat','w2');

beta1 = zeros(nSim, ns, nweek);
for iter = 1:nSim
    season1_data = result_beta{iter}{1};
    for week = 1:nweek
        beta1(iter, :, week) = season1_data(:, week);
    end
end
% save('result/beta1.mat','beta1');
beta2 = zeros(nSim, ns, nweek);
for iter = 1:nSim
    season2_data = result_beta{iter}{2};
    for week = 1:nweek
        beta2(iter, :, week) = season2_data(:, week);
    end
end
% save('result/beta2.mat','beta2');

init1 = zeros(32, nSim);
init2 = zeros(32, nSim);
for i = 1:nSim
    cell_array = mean_result_save{i};  % mean_result_save：拟合得到result_save后对第二个维度（粒子数10000）求均值后
    matrix1 = cell_array{1};
    matrix2 = cell_array{2};
    init1(:, i) = matrix1(:, 1);  % 取第一周数据作为init
    init2(:, i) = matrix2(:, 1);
end
% save('result/init1.mat','init1');
% save('result/init2.mat','init2');

init_param = stoc_params();

%% forecasting test
% 调用 Stoch_Iteration 函数
[sym_incidence1, mean_sym_incidence1, std_sym_incidence1, check_beta1] = ...
    replay(beta1, w1, init_param, init1);

[sym_incidence2, mean_sym_incidence2, std_sym_incidence2, check_beta2] = ...
    replay(beta2, w2, init_param, init2);

%% save forecasting result
mean_sym_incidence = zeros(nweek,2);
std_sym_incidence = zeros(nweek,2);
mean_sym_incidence1 = reshape(mean_sym_incidence1, [], 1);
std_sym_incidence1 = reshape(std_sym_incidence1, [], 1);
mean_sym_incidence2 = reshape(mean_sym_incidence2, [], 1);
std_sym_incidence2 = reshape(std_sym_incidence2, [], 1);
mean_sym_incidence(:,1) = mean_sym_incidence1;
mean_sym_incidence(:,2) = mean_sym_incidence2;
std_sym_incidence(:,1) = std_sym_incidence1;
std_sym_incidence(:,2) = std_sym_incidence1;
save('forecast/mean_sym_incidence.mat', 'mean_sym_incidence');
save('forecast/std_sym_incidence.mat', 'std_sym_incidence');

Observation = stateNums(:,1)'*mu;
Observation = Observation(1:nweek)';
mat_data = load('forecast/mean_values.mat','mean_values');
fit_data = mat_data.mean_values(1:nweek, 1);

figure;
hold on;
plot(fit_data);
plot(mean_sym_incidence1);
plot(Observation);
hold off;
saveas(gcf, strcat('forecast/',num2str(ns),'plot.png'));

fit_beta = mean(beta1,1);
fit_beta = mean(fit_beta,2);
fit_beta = squeeze(fit_beta);
save(strcat('forecast/',num2str(ns),'fit_beta.mat'),'fit_beta');

forecast_beta = mean(check_beta1,2);
forecast_beta = squeeze(forecast_beta);
save(strcat('forecast/',num2str(ns),'forecast_beta.mat'),'forecast_beta');
