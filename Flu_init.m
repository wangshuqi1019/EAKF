function prior = Flu_init(obs, P, param)
    % init the particle at time step 1.
    % x(1) = S, x(2) = E, x(3) = Y, x(4) = A, x (5) = R, x(6) = H, x(7) = D
    % four age groups
    % x(29) = symptomatic_incidence, x(30) = asymptomatic_incidence,
    % x(31) = x(25)+x(26), x(32) = beta
    prior = zeros(P,32);
    
    Y_mean = obs;
    
    prior(:,1) = param.init_S0(1); % S
%     prior(:,2) = 0;
%     prior(:,3) = Y_mean * param.pop(1);
%     prior(:,4) = 0;
%     prior(:,5) = param.pop(1)-prior(:,1)-prior(:,3);
    prior(:,2) = Y_mean * param.pop(1);
    prior(:,3) = param.sigma(1) * Y_mean * param.pop(1);
    prior(:,4) = (1-param.sigma(1)) * Y_mean * param.pop(1);
    prior(:,5) = param.pop(1)-prior(:,1)-prior(:,2);
    prior(:,6) = 0;
    prior(:,7) = 0;
    
    prior(:,8) = param.init_S0(2); % S
%     prior(:,9) = 0;
%     prior(:,10) = Y_mean * param.pop(2);
%     prior(:,11) = 0;
%     prior(:,12) = param.pop(2)-prior(:,8)-prior(:,10);
    prior(:,9) = Y_mean * param.pop(2);
    prior(:,10) = param.sigma(2) * Y_mean * param.pop(2);
    prior(:,11) = (1-param.sigma(2)) * Y_mean * param.pop(2);
    prior(:,12) = param.pop(2)-prior(:,8)-prior(:,9);
    prior(:,13) = 0;
    prior(:,14) = 0;

    prior(:,15) = param.init_S0(3); % S
%     prior(:,16) = 0;
%     prior(:,17) = Y_mean * param.pop(3);
%     prior(:,18) = 0;
%     prior(:,19) = param.pop(3)-prior(:,15)-prior(:,17);
    prior(:,16) = Y_mean * param.pop(3);
    prior(:,17) = param.sigma(3) * Y_mean * param.pop(3);
    prior(:,18) = (1-param.sigma(3)) * Y_mean * param.pop(3);
    prior(:,19) = param.pop(3)-prior(:,15)-prior(:,16);
    prior(:,20) = 0;
    prior(:,21) = 0;
    
    prior(:,22) = param.init_S0(4); % S
%     prior(:,23) = 0;
%     prior(:,24) = Y_mean * param.pop(4);
%     prior(:,25) = 0;
%     prior(:,26) = param.pop(4)-prior(:,22)-prior(:,24);
    prior(:,23) = Y_mean * param.pop(4);
    prior(:,24) = param.sigma(4) * Y_mean * param.pop(4);
    prior(:,25) = (1-param.sigma(4)) * Y_mean * param.pop(4);
    prior(:,26) = param.pop(4)-prior(:,22)-prior(:,23);
    prior(:,27) = 0;
    prior(:,28) = 0;
    
    prior(:,29) = param.sigma(1)*Y_mean;
    prior(:,30) = (1-param.sigma(1))*Y_mean;
    prior(:,31) = Y_mean;

    prior(:,32) = param.init_betamnl_mean + param.init_betamnl_var * rand;
end
