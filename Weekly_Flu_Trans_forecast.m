function Xn = Weekly_Flu_Trans_forecast(X, param)
    Xn = X;
    weekly_infection_count = 0;
    weekly_symptomatic_count = 0;
    weekly_asymptomatic_count = 0;
    for i = 1 : 7
        Xn = Flu_Trans_forecast(Xn, param);
        weekly_symptomatic_count = weekly_symptomatic_count + Xn(29);
        weekly_asymptomatic_count = weekly_asymptomatic_count + Xn(30);
        weekly_infection_count = weekly_infection_count + Xn(31);
    end
    Xn(29) = weekly_symptomatic_count;
    Xn(30) = weekly_asymptomatic_count;
    Xn(31) = weekly_infection_count;
end

