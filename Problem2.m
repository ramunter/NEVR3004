
figure(1);
for i  = 20:29
        subplot(5, 2, i-19); 
        plot(plottingangles,firing_rate_matrix(i,:), 'lineWidth', 2);20:30
        hold on;
        [gaussian_fit, gaus_gof] = fit(plottingangles', firing_rate_matrix(i,:)', 'gauss1');
        plot(gaussian_fit);
        hold on;
        [cos_fit, cos_gof] = fit(plottingangles', firing_rate_matrix(i,:)', 'sin1');   
        plot(cos_fit);
end

figure(2);
for i = 20:29
        [gaussian_fit, gaus_gof] = fit(plottingangles', firing_rate_matrix(i,:)', 'gauss1');
        [cos_fit, cos_gof] = fit(plottingangles', firing_rate_matrix(i,:)', 'sin1');   
        scatter(gaus_gof.rmse, cos_gof.rmse);
        hold on;
end