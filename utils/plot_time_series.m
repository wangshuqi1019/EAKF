function [] = plot_time_series(data, param)
    % plot_time_series
    % data: [Ny,T,B] 一般Ny == 1
    
%     figure;
    hold on
    plot(mean(data,3), 'r');
       
    if param.plot_quantile
        uci = param.upper_confidence_interval;
        lci = param.lower_confidence_interval;
        
        plot(quantile(data, uci, 3), 'k');
        plot(quantile(data, lci, 3), 'k');
        legend('mean', strcat(num2str(uci*100), " % interval"), strcat(num2str(lci*100), "% interval")); 
    else
        mean_y = mean(data,3);
        std_y = std(data, 0, 3);
        plot(mean_y + param.nstd * std_y, 'k');
        plot(mean_y - param.nstd * std_y, 'k');
        legend('mean', strcat("upper ", num2str(param.nstd), " std"), strcat("lower ", num2str(param.nstd), " std")); 
    end
%     xlim([42,50])
end


