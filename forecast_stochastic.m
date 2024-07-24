% version from Zengyang
function [xk_forcast,yh_forecast,beta_forecast] = forecast_stochastic(pf, param, forcast_setting, flag_plot)
% filter_start_time = replay_setting.filter_start_time;
% filter_end_time = replay_setting.filter_end_time;

obs = forcast_setting.obs;
pf_forecast = pf;
param_forecast = param;
ny = forcast_setting.ny;
nx = size(pf_forecast.particles,1); 

% S_forecast = zeros(ny, forcast_setting.T, forcast_setting.NIteration); 
% I_forecast = zeros(ny, forcast_setting.T, forcast_setting.NIteration); 
% A_forecast = zeros(ny, forcast_setting.T, forcast_setting.NIteration); 
% R_forecast = zeros(ny, forcast_setting.T, forcast_setting.NIteration); 
% H_forecast = zeros(ny, forcast_setting.T, forcast_setting.NIteration); 
% D_forecast = zeros(ny, forcast_setting.T, forcast_setting.NIteration); 
xk_forcast = zeros(11,forcast_setting.T, forcast_setting.NIteration);
yh_forecast = zeros(ny, forcast_setting.T, forcast_setting.NIteration); 
beta_forecast = zeros(ny, forcast_setting.T, forcast_setting.NIteration); 
for iter = 1 : forcast_setting.NIteration
    xkm1 = forcast_setting.init_x0;

    for k = forcast_setting.filter_start_time : forcast_setting.filter_end_time
%         disp(k)
        pf_forecast.k = k;
        param_forecast.ScalingBeta = forcast_setting.ScalingBeta;
    
        param_k = k;
        % draw one from distribution
        idx = randsample(1:pf_forecast.Ns, 1, true, pf_forecast.w(:, param_k));
        param_curr = pf_forecast.particles(:, idx, param_k);

        xkm1(10:11) = param_curr(10:11);
        xk = Weekly_Flu_Trans_forecast(param_k, xkm1, param_forecast, 0);
        
        xkm1 = xk;
        yh_forecast(:,k-forcast_setting.filter_start_time+1,iter) = obs(k, xk, 0);
        xk_forcast(:,k-forcast_setting.filter_start_time+1,iter) = xk;
        
        beta_forecast(:,k-forcast_setting.filter_start_time+1,iter) = (forcast_setting.ScalingBeta* param_curr(10)*param.delta1*(param.alpha+param.delta2))/(param.sigma*param.delta1+(1-param.sigma)*param.delta2);
    end

end


end

