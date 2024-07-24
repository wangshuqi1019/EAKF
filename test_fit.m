clear;
addpath('data_shenzhen/')

%% load data
nperiod = 2;
nweek = 56;
[num1,~,~] = xlsread('2017-2019.xlsx');
stateNums = zeros(nweek*nperiod,1); % 56 weeks, 2 years
temp1 = num1(:,5);
stateNums(1:length(temp1)) = temp1;

%% load parameters
nSim = 50;  % origin 100, repeating times, limited by server
ns = 1000;

iter = 0;
mu = 0.1;  % ILI+ mapping
init = 0.65;  % init S0
param_variance = 0.01;
betavar = 0.001;
betamean = 0.02;
oev = 1e-5;
% 直接预测beta，没有用R0计算beta
load('result/S0_tri.mat');  % based on init=0.65
% for param_variance = 0.005:0.005:0.05
%     for mu = 0.1:0.1:0.1
%         for betamean = 0.01:0.005:0.02
%             for betavar = 0.001:0.0003:0.0016
%                 for oev = 1e-6:5e-6:1e-4
%                     for init = 0.6:0.05:0.7
%                         for seed = 300:20:380
%                             rng(seed);
                            %% S0
%                             n_samples = nSim;
%                             lhs_samples = lhsdesign(n_samples, 2);
%                             s1_samples = init + lhs_samples(:, 1) * 0.1;
%                             s2_samples = init + lhs_samples(:, 2) * 0.1;
%                             S0_tri = zeros(nSim,2);
%                             for i = 1:n_samples
%                                 s1 = s1_samples(i);
%                                 s2 = s2_samples(i);
%                                 S0_tri(i,1) = s1;
%                                 S0_tri(i,2) = s2;   
%                             end
                            %% test
                            iter=iter+1;
                            Observation = stateNums(:,1)'*mu;
                            [MSE_list, result_save, result_w, result_beta] = ...
                                Flu_dif_period_test(stateNums, mu, S0_tri(1:nSim,:), param_variance, betamean, betavar, oev, ns);
                            mse = sum(MSE_list);
                            disp(mse);
%                             disp(['var = ' num2str(param_variance) ', betamean = ' num2str(betamean) ', betavar = ' num2str(betavar) ', oev = ' num2str(oev) ', mse = ' num2str(mse)]);
                            %% save results
%                             if sum(MSE_list)<0.003
                                dir_name = strcat('result/iter',num2str(iter));
                                mkdir(dir_name);
%                             save(strcat(dir_name,'/S0_tri.mat'),'S0_tri');                                
                                
                                para.param_variance=param_variance;
                                para.mu=mu;
                                para.betamean=betamean;
                                para.betavar=betavar;
                                para.oev=oev;
                                para.init=init;
%                                 save('result/para.mat','para');
                                save(strcat(dir_name,'/para.mat'),'para');
                                
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
%                                 save('result/mean_values.mat','mean_values');
%                                 save('result/std_values.mat','std_values');
                                save(strcat(dir_name,'/mean_values.mat'),'mean_values');
                                save(strcat(dir_name,'/std_values.mat'),'std_values');
%                                 save(strcat(dir_name,'/result_w.mat'),'result_w');
%                                 save(strcat(dir_name,'/result_beta.mat'),'result_beta');

                                plot_obs_fit(Observation, mean_values, std_values);
                                saveas(gcf, strcat('result/figure/',num2str(iter),'plot.png'));
%                             end
                        clear MSE_list
                        clear para
                        clear -regexp ^result
                        clear -regexp ^mean
                        clear -regexp ^std
                        clear -regexp ^all
                        clear -regexp ^data
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end

