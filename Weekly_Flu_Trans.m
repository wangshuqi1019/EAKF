function Xn = Weekly_Flu_Trans(l, X, param, vk)
    Xn = X;
    
    Xn(32) = Xn(32) + param.param_variance * randn;
    if Xn(32)<0
        Xn(32)=0.01;
    end
    
    weekly_infection_count = 0;
    weekly_symptomatic_count = 0;
    weekly_asymptomatic_count = 0;
    for i = 1 : 7
        Xn = Flu_Trans(param.week_idx * 7 + i, Xn, param, vk);
        weekly_symptomatic_count = weekly_symptomatic_count + Xn(29);
        weekly_asymptomatic_count = weekly_asymptomatic_count + Xn(30);
        weekly_infection_count = weekly_infection_count + Xn(31);
    end
    Xn(29) = weekly_symptomatic_count;
    Xn(30) = weekly_asymptomatic_count;
    Xn(31) = weekly_infection_count;
    
    if Xn(32)<0
        Xn(32)=0.01;
    end
end

