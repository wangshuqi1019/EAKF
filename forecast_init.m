% version from Zengyang
function prior = forecast_init(S, obs, param)
    % init the particle at time step 1.
    % x(1) = S, x(2) = I, x(3) = A, x(4) = R, x (5) = H, x(6) = D
    % x(7) = symptomatic_incidence, x(8) = asymptomatic_incidence,
    % x(9) = x(7)+x(8), x(10)=R0mxl, x(11) = R0mnl
    prior = zeros(1,11);
    
    Y_mean = obs;
    
    prior(:,1) = S; % S
%     disp(prior(:,1));
%     disp(prior(:,1))
    prior(:,2) = param.sigma * Y_mean;
    prior(:,3) = (1-param.sigma) * Y_mean;
    prior(:,4) = 1-prior(:,1)-prior(:,2)-prior(:,3);
    prior(:,5) = 0;
    prior(:,6) = 0;
    prior(:,7) = param.sigma * Y_mean;
    prior(:,8) = (1-param.sigma) * Y_mean;
    prior(:,9) = Y_mean;
    prior(:,10) = param.init_R0mxl_mean + param.init_R0mxl_var * rand;
    prior(:,11) = param.init_R0mnl_mean + param.init_R0mnl_var * rand;

end