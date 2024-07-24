function [pf, param] = Flu_filter(Observation, param_variance, init_S0, betamean, betavar, oev, ns)
% hyper_param.param_variance = param_variance;
% hyper_param.init_S0 = init_S0;

%% Process equation x[k] = sys(k, x[k-1], u[k]);
nx = 32;  % number of states
sys = @Weekly_Flu_Trans;
observation_idx = 29;  % symptomatic incidence

%% Observation equation y[k] = obs(k, x[k], v[k]);
ny = 1;  % number of observations
obs = @(k, xk, vk) xk(observation_idx);   % (returns column vector)

%% PDF of process noise and noise generator function
nu = 1;  % size of the vector of process noise
sigma_u = 0.0001;  % T0-BE-Determine

%% Initial PDF
gen_x0 = @Flu_init;  % sample from p_x0 (returns column vector)

%% Observation likelihood PDF p(y[k] | x[k])
p_yk_given_xk = @Likelihood_Emission;

%% Number of time steps
T = length(Observation);

%% Separate memory space
nv = 1; % dummy
gen_sys_noise = @(u) normrnd(0, sigma_u);  % dummy

% x, y, u, v 分别是状态向量、观测向量、过程噪声向量和观测噪声向量的存储矩阵
x = zeros(nx,T);  y = zeros(ny,T);
u = zeros(nu,T);  v = zeros(nv,T);

%% Epidemic Model Parameter
% % covariate data 
% ah_begin_week = 16;

% Parameter
param.sigma = 0.55*ones(1,4);
param.omega = 0.36*ones(1,4);
param.p_eya = 2/3*ones(1,4);
param.p_ar = 1/3*ones(1,4);
param.p_yr = 1/5*(1-0.0146)*ones(1,4);
param.p_hr = 0.9981*ones(1,4);
% param.p_yh = 0.0146;
param.p_yh = [0.0070, 0.0027, 0.0083, 0.0909];
param.p_rs = 1/(3*365)*ones(1,4);
% param.p_hd = 0.0019;
param.p_hd = [0.000050, 0.000072, 0.000595, 0.001570];

param.N = 17681600;  % Shenzhen population
param.pop = [0.0566, 0.1314, 0.7585, 0.0535];  % population proportions in four age groups
param.phi = [[1.409660734, 1.768548709, 3.751113208, 0.350060933];...
                        [0.346000233, 10.64471027, 6.08382839, 0.286702683];...
                        [0.570072983, 4.073960628, 10.62071562, 0.482763604];...
                        [0.396853292, 1.754034964, 3.745765744, 1.173235326]];  % age-speciofic contact matrix
                    
param.betamnl = 0; param.betamnu = 1.0;
% param.init_S0 = hyper_param.init_S0;
param.init_S0 = init_S0;

param.init_betamnl_mean = betamean;
param.init_betamnl_var = betavar;

param.obs_idx = observation_idx;
% param.param_variance = hyper_param.param_variance; % 0.3; %state 用 0.18, us 用0.3
param.param_variance = param_variance; % 0.3; %state 用 0.18, us 用0.3

%% Parameter of the algorithm: \beta = a*F_t+b (a is learned )
% Comment either one to deterine the function param.betaf = @(X,param)X(7)*param.Fatigue_PM3n_k + param.b;
% learn function param.betaf = @(X,param)X(7); % learn it freely
pf.Ns              = ns;               % number of particles
pf.k               = 1;                   % initial iteration number
pf.w               = zeros(pf.Ns, T);     % weights
pf.particles       = zeros(nx, pf.Ns, T); % particles
pf.gen_x0          = gen_x0;              % function for sampling from initial pdf p_x0
pf.p_yk_given_xk   = p_yk_given_xk;       % function of the observation likelihood PDF p(y[k] | x[k])
pf.gen_sys_noise   = gen_sys_noise;       % function for generating system noise # no used.
pf.obs = obs;

%% Estimate state
filter_start_time = 1; % filter_end_time should > 1
filter_end_time = length(Observation);

observation_curr = Observation(filter_start_time : filter_end_time);

% init
T = filter_end_time - filter_start_time + 1;
xh0 = gen_x0(observation_curr(1), 1, param);

xh = zeros(nx, T); xh(:, 1) = xh0;
yh = zeros(ny, T); yh(:, 1) = obs(1, xh0, 0);
yhs = zeros(pf.Ns, T); yhs(:, 1) = yh(:, 1);

% B = zeros(nx, T*7);

week_idx = filter_start_time - 1;

for k = 2 : filter_end_time - filter_start_time + 1
    pf.k = k;
    
%     if mod(k, 10) == 0
% %         disp(k + filter_start_time)
%     end

    oev5_factor = oev; % 可调，origin 0.2
%     oev5_factor = oev + 0.002 * randn(); % 可调，origin 0.2
    param.oev = oev5_factor * sqrt(mean(observation_curr(:, max(k-1,1):k-1).^2));
%     param.oev = param.oev + 0.0002 * randn; % 增加随机扰动
    param.oev = param.oev + 0.01 * randn; % 增加随机扰动
%     disp(param.oev)
%     param.oev = oev5_factor * mean(observation_curr(:, max(k-3,1):k-1));
    if param.oev < 1e-5 % 可调
        param.oev = 1e-5;
    end
    param.oev = 1e-5 + 1e-6*randn;
    param.week_idx = week_idx;

%     smoothed_observation = movmean(observation_curr, [2 0], 2);
    smoothed_observation = movmean(observation_curr, [1 1], 2);
    [xh(:,k), pf] = particle_filter(sys, smoothed_observation(:, k), pf, 'systematic_resampling', param); % 更新状态估计
%     [xh(:,k), pf] = particle_filter(sys, observation_curr(:, k), pf, 'systematic_resampling', param); % 更新状态估计
    % filtered observation
    yh(:,k) = obs(k, xh(:,k), 0); % 更新观测估计
    week_idx=  week_idx + 1;
end
end
