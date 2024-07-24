% version from Zengyang
function replay_mean(pf, param, replay_setting)

filter_start_time = replay_setting.filter_start_time;
filter_end_time = replay_setting.filter_end_time;
week_idx = replay_setting.week_idx;
obs = replay_setting.obs;
pf_replay = pf;
param_replay = param;
ny = replay_setting.ny;
nx = size(pf_replay.particles,1); 

yh_replay = zeros(ny, replay_setting.T); 
xkm1 = zeros(nx,1);

for i = 1:pf_replay.Ns
     xkm1 = xkm1 + pf_replay.w(i, 2) * pf_replay.particles(:, i, 2);
    % pare_idx = pf_replay.parenets(i,3);
    % xkm1 = xkm1 +  pf_replay.w(i, 3) * pf_replay.particles(:, pare_idx, 2);
end

for k = 2 : replay_setting.T-1
    pf_replay.k = k;
    param.week_idx = week_idx;    
    param_replay.ScalingBeta = replay_setting.ScalingBeta;

    param_curr = zeros(nx,1);
    param_k = k+1;
    for i = 1:pf_replay.Ns
        param_curr = param_curr + pf_replay.w(i, param_k) ...
                    * pf_replay.particles(:, i, param_k);
    end
    xkm1(3:4) = param_curr(3:4);
    xk = Weekly_Flu_Trans(k, xkm1, param_replay, 0);
    
    xkm1 = xk;
    week_idx=  week_idx + 1;

    yh_replay(:,k+1) = obs(k, xk, 0);
end


figure;
hold on
plot((filter_start_time: filter_end_time ), squeeze(yh_replay), 'r')
plot((filter_start_time: filter_end_time ), replay_setting.observation_curr, 'g');
legend('replay', 'ground-truth');
xlim([filter_start_time+1, filter_end_time]);
end
