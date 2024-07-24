function [MSE_list, result_save, result_w, result_beta] = ...
            Flu_dif_period2(stateNums, mu, S0_tri, param_variance, betamean, betavar, oev)
    addpath('utils/')

    % ILI+
    Observation = stateNums(:,1)'*mu;
    nperiod = 2; % 两年周期，每个周期内按照56周拟合
    
    MSE_list = zeros(50,1);
    %% filter transition paramete
    parpool('local', 25);
    parfor i_iter = 1:50
%         rng((i_iter+50)*10);
        start_time  = 1;
        pfs = {};
        params = {};
        save_state = {};
        MSE = 0;
        beta = {};
        w = {};
        for i = 1 : nperiod
            end_time = start_time + 55;
            init_S0 = S0_tri(i_iter,i)*[0.0566, 0.1314, 0.7584, 0.0535];  % 乘以四个年龄组的人口比例
            [pf, param] = Flu_filter(Observation(start_time:end_time), param_variance, init_S0, betamean, betavar, oev);
            pf.observation_curr = Observation(start_time:end_time);
            pfs(i) = {pf};
            save_state(i) = {pf.particles};
            % 补充保存权值
            w(i) = {pf.w};
            % 补充保存传染率
            beta(i) = {squeeze(pf.particles(32, :, :))};
            temp_data = permute(save_state{i}(29,:,:), [1,3,2]); % 观察有症状incidence
            MSE = MSE + sum((mean(temp_data,3)-Observation(start_time:end_time)).^2);
            params(i) = {param};
            start_time = start_time + 56;
        end
        MSE_list(i_iter,1) = MSE; % 衡量模型预测值与实际观测值之间的偏差
        result_save{i_iter} = save_state;
        result_w{i_iter} = w;
        result_beta{i_iter} = beta;
    end
    delete(gcp('nocreate'))   
end