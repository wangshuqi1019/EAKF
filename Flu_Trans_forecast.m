function Xn = Flu_Trans_forecast(X, param)

    pop = param.pop;
    phi = param.phi;
    
    sigma = param.sigma;
    omega = param.omega;
    p_ar = param.p_ar;
    p_yr = param.p_yr;
    p_hr = param.p_hr;
    p_yh = param.p_yh;
    p_hd = param.p_hd;
    p_rs = param.p_rs;
    p_eya = param.p_eya;
    
    Xn = X;
    beta = Xn(32);
    Xn_next = Xn;
    
    %%
    new1 = beta*Xn(1)*((Xn(3)+omega(1)*Xn(4))*phi(1,1)/pop(1)+...
                        (Xn(10)+omega(1)*Xn(11))*phi(1,2)/pop(2)+...
                        (Xn(17)+omega(1)*Xn(18))*phi(1,3)/pop(3)+...
                        (Xn(24)+omega(1)*Xn(25))*phi(1,4)/pop(4));
    Xn_next(1) = Xn(1) - new1 + p_rs(1)*Xn(5); % S_1 (age 0-5)
    Xn_next(2) = Xn(2) + new1 - p_eya(1)*Xn(2); % E_1 (age 0-5)
    Xn_next(3) = Xn(3) + sigma(1)*p_eya(1)*Xn(2) - p_yh(1)*Xn(3) - p_yr(1)*Xn(3); % Y_1 (age 0-5)
    Xn_next(4) = Xn(4) + (1-sigma(1))*p_eya(1)*Xn(2) - p_ar(1)*Xn(4); % A_1 (age 0-5)
    Xn_next(5) = Xn(5) + p_ar(1)*Xn(4) + p_yr(1)*Xn(3) + p_hr(1)*Xn(6) - p_rs(1)*Xn(5); % R_1 (age 0-5)
    Xn_next(6) = Xn(6) + p_yh(1)*Xn(3) - p_hr(1)*Xn(6) - p_hd(1)*Xn(6); % H_1 (age 0-5)
    Xn_next(7) = Xn(7) + p_hd(1)*Xn(6); % D_1 (age 0-5)
    
    new2 = beta*Xn(8)*((Xn(3)+omega(2)*Xn(4))*phi(2,1)/pop(1)+...
                        (Xn(10)+omega(2)*Xn(11))*phi(2,2)/pop(2)+...
                        (Xn(17)+omega(2)*Xn(18))*phi(2,3)/pop(3)+...
                        (Xn(24)+omega(2)*Xn(25))*phi(2,4)/pop(4));
    Xn_next(8) = Xn(8) - new2 + p_rs(2)*Xn(12); % S_2 (age 6-18)
    Xn_next(9) = Xn(9) + new2 - p_eya(2)*Xn(9); % E_2 (age 6-18)
    Xn_next(10) = Xn(10) + sigma(2)*p_eya(2)*Xn(9) - p_yh(2)*Xn(10) - p_yr(2)*Xn(10); % Y_2 (age 6-18)
    Xn_next(11) = Xn(11) + (1-sigma(2))*p_eya(2)*Xn(9) - p_ar(2)*Xn(11); % A_2 (age 6-18)
    Xn_next(12) = Xn(12) + p_ar(2)*Xn(11) + p_yr(2)*Xn(10) + p_hr(2)*Xn(13) - p_rs(2)*Xn(12); % R_2 (age 6-18)
    Xn_next(13) = Xn(13) + p_yh(2)*Xn(10) - p_hr(2)*Xn(13) - p_hd(2)*Xn(13); % H_2 (age 6-18)
    Xn_next(14) = Xn(14) + p_hd(2)*Xn(13); % D_2 (age 6-18)

    new3 = beta*Xn(15)*((Xn(3)+omega(3)*Xn(4))*phi(3,1)/pop(1)+...
                        (Xn(10)+omega(3)*Xn(11))*phi(3,2)/pop(2)+...
                        (Xn(17)+omega(3)*Xn(18))*phi(3,3)/pop(3)+...
                        (Xn(24)+omega(3)*Xn(25))*phi(3,4)/pop(4));
    Xn_next(15) = Xn(15) - new3 + p_rs(3)*Xn(19); % S_3 (age 19-59)
    Xn_next(16) = Xn(16) + new3 - p_eya(3)*Xn(16); % E_3 (age 19-59)
    Xn_next(17) = Xn(17) + sigma(3)*p_eya(3)*Xn(16) - p_yh(3)*Xn(17) - p_yr(3)*Xn(17); % Y_3 (age 19-59)
    Xn_next(18) = Xn(18) + (1-sigma(3))*p_eya(3)*Xn(16) - p_ar(3)*Xn(18); % A_3 (age 19-59)
    Xn_next(19) = Xn(19) + p_ar(3)*Xn(18) + p_yr(3)*Xn(17) + p_hr(3)*Xn(20) - p_rs(3)*Xn(19); % R_3 (age 19-59)
    Xn_next(20) = Xn(20) + p_yh(3)*Xn(17) - p_hr(3)*Xn(20) - p_hd(3)*Xn(20); % H_3 (age 19-59)
    Xn_next(21) = Xn(21) + p_hd(3)*Xn(20); % D_3 (age 19-59)

    new4 = beta*Xn(22)*((Xn(3)+omega(4)*Xn(4))*phi(4,1)/pop(1)+...
                        (Xn(10)+omega(4)*Xn(11))*phi(4,2)/pop(2)+...
                        (Xn(17)+omega(4)*Xn(18))*phi(4,3)/pop(3)+...
                        (Xn(24)+omega(4)*Xn(25))*phi(4,4)/pop(4));
    Xn_next(22) = Xn(22) - new4 + p_rs(4)*Xn(26); % S_4 (age ≥60)
    Xn_next(23) = Xn(23) + new4 - p_eya(4)*Xn(23); % E_4 (age ≥60)
    Xn_next(24) = Xn(24) + sigma(4)*p_eya(4)*Xn(23) - p_yh(4)*Xn(24) - p_yr(4)*Xn(24); % Y_4 (age ≥60)
    Xn_next(25) = Xn(25) + (1-sigma(4))*p_eya(4)*Xn(23) - p_ar(4)*Xn(25); % A_4 (age ≥60)
    Xn_next(26) = Xn(26) + p_ar(4)*Xn(25) + p_yr(4)*Xn(24) + p_hr(4)*Xn(27) - p_rs(4)*Xn(26); % R_4 (age ≥60)
    Xn_next(27) = Xn(27) + p_yh(4)*Xn(24) - p_hr(4)*Xn(27) - p_hd(4)*Xn(27); % H_4 (age ≥60)
    Xn_next(28) = Xn(28) + p_hd(4)*Xn(27); % D_4 (age ≥60)
    
    Xn_next(29) = sigma(1)*new1 + sigma(2)*new2 + sigma(3)*new3 + sigma(4)*new4; % symptomatic_incidence                                  
    Xn_next(30) = (1-sigma(1))*new1 + (1-sigma(2))*new2 + (1-sigma(3))*new3 + (1-sigma(4))*new4; % asymptomatic_incidence
    Xn_next(31) = new1+new2+new3+new4; % incidence
    
    Xn = Xn_next;
end
