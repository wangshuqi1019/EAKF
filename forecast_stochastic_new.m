% version from Zengyang
function [xk_forcast,yh_forecast,beta_forecast] = forecast_stochastic_new(pf, param, forcast_setting, flag_plot)
% filter_start_time = replay_setting.filter_start_time;
% filter_end_time = replay_setting.filter_end_time;

obs = forcast_setting.obs;
pf_forecast = pf;
param_forecast = param;
ny = forcast_setting.ny;
nx = size(pf_forecast.particles,1); 


xk_forcast = zeros(11,forcast_setting.T, forcast_setting.NIteration);
yh_forecast = zeros(ny, forcast_setting.T, forcast_setting.NIteration); 
beta_forecast = zeros(ny, forcast_setting.T, forcast_setting.NIteration); 
for iter = 1 : forcast_setting.NIteration
    xkm1 = forcast_setting.init_x0;

    for k = forcast_setting.filter_start_time : forcast_setting.filter_end_time
        pf_forecast.k = k;
        param_forecast.ScalingBeta = forcast_setting.ScalingBeta;
        if k>26     % north 26, south 48
            param_k = k-9;
        else
            param_k = k;
        end
        % draw one from distribution
        idx = randsample(1:pf_forecast.Ns, 1, true, pf_forecast.w(:, param_k));
        param_curr = pf_forecast.particles(:, idx, param_k);

        xkm1(10:11) = param_curr(10:11);
        xk = Weekly_Flu_Trans(k, xkm1, param_forecast, 0);
        
        xkm1 = xk;
        yh_forecast(:,k-forcast_setting.filter_start_time+1,iter) = obs(k, xk, 0);
        xk_forcast(:,k-forcast_setting.filter_start_time+1,iter) = xk;
        
        beta_forecast(:,k-forcast_setting.filter_start_time+1,iter) = (forcast_setting.ScalingBeta* param_curr(10)*param.delta1*(param.alpha+param.delta2))/(param.sigma*param.delta1+(1-param.sigma)*param.delta2);
    end

end

if(flag_plot)
    figure;
    hold on
    plot(mean(yh_forecast,3), 'r')
    plot(quantile(yh_forecast, 0.95, 3), 'k');
    plot(quantile(yh_forecast, 0.05, 3), 'k');
%     xlim([3, forcast_setting.T])
    legend('forecast', '90% confident level','10% confident level');
end
end

