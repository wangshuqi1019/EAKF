function param = stoc_params()
    param.N = 17681600;

    % 各年龄组人口比例
    param.pop = [0.0566; 0.1314; 0.7585; 0.0535];

    % 接触矩阵
    param.phi = [
        1.409660734, 1.768548709, 3.751113208, 0.350060933;
        0.346000233, 10.64471027, 6.08382839, 0.286702683;
        0.570072983, 4.073960628, 10.62071562, 0.482763604;
        0.396853292, 1.754034964, 3.745765744, 1.173235326
    ];

    param.age_n = 4;  % 年龄组

    % 有症状的比例
    param.sigma = [0.55; 0.55; 0.55; 0.55];

    % 无症状的相对传染性
    param.omega = [0.36; 0.36; 0.36; 0.36];

    % 离开暴露状态的概率
    param.p_eya = [1/1.5; 1/1.5; 1/1.5; 1/1.5];

    % 无症状恢复率
    param.p_ar = [1/3; 1/3; 1/3; 1/3];

    % 有症状恢复率
    param.p_yr = [1/5*(1-0.0146); 1/5*(1-0.0146); 1/5*(1-0.0146); 1/5*(1-0.0146)];

    % 住院患者恢复率
    param.p_hr = [0.9981; 0.9981; 0.9981; 0.9981];

    % 住院患者死亡率
    param.p_hd = [0.00005; 0.000072; 0.000595; 0.00157];

    % 有症状住院率
    param.p_yh = [0.007; 0.0027; 0.0083; 0.0909];

    % 离开恢复状态的速率
    param.p_rs = [1/(3*365); 1/(3*365); 1/(3*365); 1/(3*365)];

    % 疫苗有效性
    param.v_e = [1/2; 9/10; 7/10; 1/2; 3/10; 1/10];

    % 疫苗子仓室的转换率
    param.f = [1/30; 1/30; 1/30; 1/30; 1/30];

    param.param_variance = 0.005;
end
