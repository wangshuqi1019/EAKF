% version from Zengyang
function lines = replay_stochastic(pf, param, replay_setting)

filter_start_time = replay_setting.filter_start_time;
filter_end_time = replay_setting.filter_end_time;
week_idx = replay_setting.week_idx;
obs = replay_setting.obs;
pf_replay = pf;
param_replay = param;
ny = replay_setting.ny;
nx = size(pf_replay.particles,1); 
param_replay.trans_variance = 0.01;

yh_replay = zeros(ny, replay_setting.T, replay_setting.NIteration);
beta_replay = zeros(ny, replay_setting.T, replay_setting.NIteration);
pf_replay.w(:, replay_setting.T+1) = 1/pf_replay.Ns;
for iter = 1 : replay_setting.NIteration
    xkm1 = zeros(nx,1);
    for i = 1:pf_replay.Ns
        xkm1 = xkm1 + pf_replay.w(i, 2) * pf_replay.particles(:, i, 2);
        % pare_idx = pf_replay.parenets(i,3);
        % xkm1 = xkm1 +  pf_replay.w(i, 3) * pf_replay.particles(:, pare_idx, 2);
    end
    yh_replay(:, 2, iter) = obs(0, xkm1, 0);
    beta_replay (:, 2, iter) = xkm1(4)/param.D;

    for k = 2: replay_setting.T-1

        pf_replay.k = k;
        param.week_idx = week_idx;    

%         if k>=5
%             param_replay.ScalingBeta = replay_setting.ScalingBeta*1.2;
%         else
%             param_replay.ScalingBeta = replay_setting.ScalingBeta;
%         end
        
   
        param_k = k;
        
        % draw one from distribution
        idx = randsample(1:pf.Ns, 1, true, pf_replay.w(:, param_k));
        param_curr = pf_replay.particles(:, idx, param_k);


        xkm1(3:4) = param_curr(3:4);
        xk = Weekly_Flu_Trans(k, xkm1, param_replay, 0);
        
        xkm1 = xk;
        week_idx=  week_idx + 1;
    
        yh_replay(:, k, iter) = obs(k, xk, 0);
        beta_replay (:, k, iter) = xkm1(4)/param.D;
    end
end

% figure;
hold on
plot((filter_start_time: filter_end_time ), mean(yh_replay,3), 'r');
plot((filter_start_time: filter_end_time ), replay_setting.observation_curr, 'g');

plot((filter_start_time: filter_end_time ), quantile(yh_replay, 0.95, 3), 'k');
plot((filter_start_time: filter_end_time ), quantile(yh_replay, 0.05, 3), 'k');
legend('replay', 'ground-truth', "95% interval", "5% interval");
xlim([filter_start_time+1, filter_end_time]);
lines = {yh_replay, replay_setting.observation_curr,beta_replay};

end

