function [sym_incidence, mean_sym_incidence, std_sym_incidence, check_beta] = ...
            replay(init_beta, init_w, param, init)
    %%
    result = cell(1,50);
    check_beta = zeros(56,50);
    
    xk_forecast = zeros(32, 56, 50);  % 57个状态，52周，100次模拟
    
    % 用struct存储每次实验的结果
    for i_iter = 1:50
        % 已知迭代次数*粒子数*周期数大小的二维数组存储的传染率result_beta和权值result_w
        pf_forecast.beta = squeeze(init_beta(i_iter, :, :));
        pf_forecast.w = squeeze(init_w(i_iter, :, :));
%         all_beta = (pf_forecast.beta).*(pf_forecast.w);
%         all_beta = squeeze(sum(all_beta,1));
        % 已知初始化时刻各状态数据init
        xkm1 = init(:, i_iter);
%         xk_forecast(:, 1, i_iter) = xkm1;
       
        for k = 1:52
            param_k = k;
%             beta = all_beta(param_k);
            weights = pf_forecast.w(:, param_k+1);
            beta_values = pf_forecast.beta(:, param_k+1);
%             weights = pf_forecast.w(:, param_k);
%             beta_values = pf_forecast.beta(:, param_k);
            beta = sum(weights .* beta_values);
%             weights = pf_forecast.w(:, param_k);
            xkm1(32) = beta;
            % 开始进行迭代，包括所有粒子所有天
            xk = Weekly_Flu_Trans_forecast(xkm1, param);
            xkm1 = xk;
            xk_forecast(:, k, i_iter) = xk;
            check_beta(k,i_iter)=beta;
        end
        
        result{i_iter} = pf_forecast;
    end
    
    sym_incidence = squeeze(xk_forecast(29, :, :));
    mean_sym_incidence = mean(sym_incidence, 2);
    std_sym_incidence = std(sym_incidence, 0, 2);
    sym_incidence = sum(mean_sym_incidence);
end
