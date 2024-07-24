function draw_percentiles(data, list_pcr, x_axis)
    for i = 1 : length(list_pcr)
        plot(squeeze(x_axis), squeeze(prctile(data, list_pcr(i), 1)));
    end
end

