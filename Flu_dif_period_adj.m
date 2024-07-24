function [flag] = Flu_dif_period_adj(statePop,stateNums,states_26_names,period_len,real_stateNums)
    addpath('utils/')
%     load stateInfos_real_autumn.mat
    dir_name = strcat('result/result_',int2str(period_len),'_different_S');
    mkdir(dir_name);
    parpool(8)
    state_index = 1;
    % ILI+
    Observation = stateNums(:,state_index)'*0.3;
%     S0_mean = [0.84,0.82,0.86,0.74,0.85]; %temp setting for weeks 36
    S0_mean = [0.84,0.82,0.84];
    S0_lower = S0_mean-0.1;
    S0_upper = S0_mean+0.1;
    nperiod = 2;
    S0_tri = zeros(100,2);
    for i = 1:nperiod
%         S0_tri(:,i) = triangular_S0(S0_mean(i),S0_lower(i),S0_upper(i));
        S0_tri(:,i) = S0_mean(i);
    end
%     S0_tri(:,:) = 0.95;
    param_variance = 1.2;
    %% filter transition paramete  
    MSE_list = zeros(100,1);
    parfor i_iter = 1:100
        start_time  = 1;
        pfs = {};
        params = {};
        save_state = {};
        MSE = 0;
        for i = 1 : nperiod
            end_time = start_time + 55;
            init_S0 = S0_tri(i_iter,i);
            [pf, param] = Flu_filter(Observation(start_time:end_time), param_variance,init_S0);
            pf.observation_curr = Observation(start_time:end_time);
            pfs(i) = {pf};
            save_state(i) = {pf.particles};
            temp_data = permute(save_state{i}(7,:,:), [1,3,2]);
            MSE = MSE + sum((mean(temp_data,3)-Observation(start_time:end_time)).^2);
            params(i) = {param};
            start_time = start_time + 56;
        end
        MSE_list(i_iter,1) = MSE;
        result{i_iter} = pfs;
        result_save{i_iter} = save_state;
        result_params{i_iter} = params;  
    end
    delete(gcp('nocreate'))
    
%     mean_forecast = {};
    fit_result = zeros(2,11,56,100);
    for i_iter = 1:100
        save_stage = result_save{i_iter};
        for k = 1:nperiod
            temp_save_stage = save_stage{k};
            for stage_index = 1:11
                mean_save_stage = permute(temp_save_stage(stage_index,:,:),[1,3,2]);
                mean_save_stage = mean(mean_save_stage,3);
                fit_result(k,stage_index,:,i_iter) = mean_save_stage;
            end
        end
    end
%     save(strcat(dir_name,'/MES_fitting_result.mat'),'MSE_list')
%     save(strcat(dir_name,'/mean_autumn_fitting_result.mat'),'mean_forecast')
%     save(strcat(dir_name,'/autumn_fitting_result_params.mat'),'result_params')
    save (strcat(dir_name,'/fit_result.mat'),'fit_result');
    mean_forecast = {};
    for i_iter = 1:100
        save_state = result_save{i_iter};
        if i_iter==1
            for k = 1:2
    %             mean_forecast{k} = {(cell2mat(mean_forecast{k})+cell2mat(save_state{k}))/2};
                mean_forecast{k} = save_state{k}/100;
            end
        else
            for k = 1:2
    %             mean_forecast{k} = {(cell2mat(mean_forecast{k})+cell2mat(save_state{k}))/2};
                mean_forecast{k} = mean_forecast{k}+save_state{k}/100;
            end
        end
    end
%     save(strcat(dir_name,'/MES_fitting_result.mat'),'MSE_list')
    save(strcat(dir_name,'/mean_autumn_fitting_result.mat'),'mean_forecast')
    save(strcat(dir_name,'/autumn_fitting_result_params.mat'),'result_params')

    start_time  = 1;
    for i = 1:nperiod
        end_time = start_time + 55;
        plot_param.plot_quantile = false;
        plot_param.upper_confidence_interval = 0.975;
        plot_param.lower_confidence_interval = 0.025;
        plot_param.nstd = 1.96;
        subplot(2,1,i);
        y_forecast = permute(mean_forecast{i}(7,:,:), [1,3,2]);
    %     subplot(2,2,i);
        for week_num = 1:55
            y_forecast(:,week_num,:) = y_forecast(:,week_num+1,:);
        end
        plot_time_series(y_forecast, plot_param);
%         title(file_name);
        xlabel('week idx');
        ylabel('incidence');
        plot(Observation(start_time:end_time), 'g');
    %     plot(Observation(start_time:end_time), 'g');
        start_time = start_time + 56;
    end
    
%     saveas(gcf,strcat(dir_name,'/fitting.png'));
    saveas(gcf,strcat(dir_name,'/fitting.fig'));

 %% forecast
    for i_iter=1:100
      pfs = result{i_iter};
      params = result_params{i_iter};
      pf_mean = pfs{1};
      param = params{1};
      for i = 2 : nperiod
          pf_mean.Ns = pf_mean.Ns+ pfs{i}.Ns;
          pf_mean.particles = cat(2, pf_mean.particles, pfs{i}.particles);
          pf_mean.w = cat(1, pf_mean.w, pfs{i}.w);
      end
      weight_2022 = 1;
      for index_2022 = 1:period_len
          pf_mean.w(1:20000,index_2022) = (1-weight_2022)*pf_mean.w(1:20000,index_2022)/4;
          pf_mean.w(20001:30000,index_2022) = weight_2022*pf_mean.w(20001:30000,index_2022);
      end
      for index_2022 = period_len+1:56
          pf_mean.w(1:20000,index_2022) = pf_mean.w(1:20000,index_2022)/4;          
          pf_mean.w(20001:30000,index_2022) = 0;
      end
      if i_iter==1
          pf = pf_mean;
      else
          pf.particles = (pf.particles+pf_mean.particles)/2;
          pf.w = (pf.w + pf_mean.w)/2;
      end
    end

    forcast_setting.obs = pf.obs;
    forcast_setting.ny = 1;
    forcast_setting.week_idx = 0;
    forcast_setting.T = 51; 
    forcast_setting.obs = pf.obs;
    forcast_setting.filter_start_time = 5; 
    forcast_setting.filter_end_time = 56;
    forcast_setting.NIteration = 200;
    forcast_setting.ScalingBeta = 1; % set the forecast ScalingBeta here!
    
    
%     potential_S = linspace(0.6,0.9,4);
    potential_S = 0.84;
    prtential_beat = 1;
    list_yh_forecast = {};
    list_beta_forecast = {};
    list_xh_forecast = {};
    list_yh_forecast_S0 = {};
    list_beta_forecast_S0 ={};
    list_xh_forecast_S0 = {};
    
    R_squared = zeros(length(potential_S),1);
    MSE = zeros(length(potential_S),1);
    
    % determain the S0
    for i = 1 : length(potential_S)
        % init S, I, R, R0mx, R0min, Obs (latter 4 is not necessary)
        if period_len ==0
           forcast_setting.init_x0 = forecast_init(potential_S(i),0.0001, param); 
        else
            forcast_setting.init_x0 = forecast_init(potential_S(i),real_stateNums(6), param);
        end
        % [yh_forecast] = forecast_mean(pf, param, forcast_setting);
        flag_plot = false;
        [xk_forcast,yh_forecast,beta_forecast] = forecast_stochastic(pf, param, forcast_setting, flag_plot);
%         MSE(i,1) = sum((Observation(208+1:208+period_len,:)'-mean(yh_forecast(1,208+1:208+period_len,:),3)).^2)/length(Observation(208+1:208+period_len,:));
        list_xh_forecast_S0(i) = {xk_forcast};
        list_yh_forecast_S0(i) = {yh_forecast};
        list_beta_forecast_S0(i) = {beta_forecast};
        if (flag_plot)
            title(strcat(states_26_names{state_index}, 'forecast with S:', num2str(potential_S(i))));
        end
    end
    
    %
%     save(strcat(dir_name,'/Texas_MSE.mat'),'MSE')
    file_name = strcat(dir_name,'/',states_26_names{state_index}, ' forecast_autumn_xh');
    save(file_name, 'list_xh_forecast_S0');
    file_name = strcat(dir_name,'/',states_26_names{state_index}, ' forecast_autumn');
    save(file_name, 'list_yh_forecast_S0');
    file_name = strcat(dir_name,'/',states_26_names{state_index}, ' forecast_autumn_beta');
    save(file_name, 'list_beta_forecast_S0');

    figure;
    hold on;
    legends = {};
    
    for i = 1 : length(potential_S)
        yh_forecast = list_yh_forecast_S0{i};
        
        if length(size(yh_forecast)) == 2
            plot(squeeze(yh_forecast))
        else
            plot(mean(yh_forecast,3))
        end
        legends(i) = {strcat(states_26_names{state_index}, ': forecast with S size', num2str(potential_S(i)))};
    end
      
    plot(real_stateNums(5:length(real_stateNums)),'LineWidth',3)
    legends(length(potential_S)+1) = {'real'};
    legend(legends,'Location','northwest')
%     xlim([3, 27])
    saveas(gcf,strcat(dir_name,'/forecast with different S.png'));
    saveas(gcf,strcat(dir_name,'/forecast with different S.fig'));

    flag = 1;
end

