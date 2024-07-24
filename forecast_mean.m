% version from Zengyang
function [yh_forecast] = forecast_mean(pf, param, forcast_setting, flag_plot)
% filter_start_time = replay_setting.filter_start_time;
% filter_end_time = replay_setting.filter_end_time;

obs = forcast_setting.obs;
pf_forecast = pf;
param_forecast = param;
ny = forcast_setting.ny;
nx = size(pf_forecast.particles,1); 

yh_forecast = zeros(ny, forcast_setting.T); 

% xkm1 = zeros(nx,1);
% for i = 1:pf_forecast.Ns
%      xkm1 = xkm1 + pf_replay.w(i, 2) * pf_replay.particles(:, i, 2);
%     % pare_idx = pf_replay.parenets(i,3);
%     % xkm1 = xkm1 +  pf_replay.w(i, 3) * pf_replay.particles(:, pare_idx, 2);
% end

xkm1 = forcast_setting.init_x0;

for k = 2 : forcast_setting.T-1
    pf_forecast.k = k;
    param_forecast.ScalingBeta = forcast_setting.ScalingBeta;

    param_curr = zeros(nx,1);
    param_k = k+1;
    for i = 1:pf_forecast.Ns
        param_curr = param_curr + pf_forecast.w(i, param_k) ...
                    * pf_forecast.particles(:, i, param_k);
    end
    xkm1(3:4) = param_curr(3:4);
    xk = Weekly_Flu_Trans(k, xkm1, param_forecast, 0);
    
    xkm1 = xk;
    yh_forecast(:,k+1) = obs(k, xk, 0);
end

if (flag_plot)
    figure;
    hold on
    plot(squeeze(yh_forecast), 'r')
    legend('forecast');
    xlim([3, forcast_setting.T])
end
end

