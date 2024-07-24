% version from Zengyang
period_length = linspace(0,47,48);
% period_length = linspace(47,47,1);
addpath('data\')
[num2,~,~] = xlsread('north_from_14_2022_2023.csv');
real_incidence = num2(:,2)*0.001;


for index = 1:length(period_length)
    
    data_address = strcat('result-2-28/result_',int2str(period_length(index)));
    load(strcat(data_address,'/mean_autumn_fitting_result.mat'))
    fit_result = mean_forecast{5};
    incidence_fit_all = permute(fit_result(9,:,:), [1,3,2]);
    Rt_fit_all = permute(fit_result(11,:,:),[1,3,2]);
    period = repmat(period_length(index),53,1);
    incidence_fit = Cal_mean_CI(incidence_fit_all);
    Rt_fit = Cal_mean_CI(Rt_fit_all);
    incidence_fit = [period,incidence_fit];
    Rt_fit = [period,Rt_fit];
    for i =1:52
        incidence_fit(i,:) = incidence_fit(i+1,:);
        Rt_fit(i,:) = Rt_fit(i+1,:);
    end
    if index==1
        save_incidence = incidence_fit;
        save_Rt = Rt_fit;
    else
        save_incidence = [save_incidence;incidence_fit];
        save_Rt = [save_Rt;Rt_fit];
    end
    save_incidence(save_incidence<0)=0;
    save_Rt(save_Rt<0)=0;
end
csvwrite('north_fit_incidence.csv',save_incidence);
csvwrite('north_fit_Rt.csv',save_Rt);
% hold on;
% plot(incidence_fit(:,1),'r');
% plot(incidence_fit(:,2), 'k');
% plot(incidence_fit(:,3), 'k');
% plot(real_incidence,'g');