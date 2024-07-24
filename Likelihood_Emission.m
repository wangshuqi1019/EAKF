% version from Zengyang
% P(X|Y)
function p = Likelihood_Emission(k, Y, X, param)

    if any(X(1:2)<0)
        p = 0;
    
    else
        error = Y - X(param.obs_idx);
        p = normpdf(error, 0, param.oev);      
    end


%     if X(29) < param.R0mnl || X(29) > param.R0mnu
%         p = 0;
%     end

    if p < 1e-20
        p = 1e-20;
    end
end