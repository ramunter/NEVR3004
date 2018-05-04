%[cellnames, spiketimes, binnedAwakeAngleData,resclu] = Import_data("Mouse12-120806/Mouse12-120806",13);

figure
xlabel('Head angle from 0-360 degrees')
ylabel('Normalized firing rate')
j=1;
num_cols = 10;
num_rows = ceil(length(cellnames)/num_cols);
firing_rate_matrix = zeros(length(cellnames),40);
for i = 1:length(cellnames)
        [firingRate,plottingangles,anglesAtFiring] = Calculate_Firing_Rate(spiketimes(cellnames{i}),binnedAwakeAngleData);
        firing_rate_matrix(i,:) = firingRate;
        subplot(num_rows, num_cols, i);
        
        normalized_firing_rate = normr(firing_rate_matrix(i,:));
        if sum(normalized_firing_rate>0.3) ~= 0
                plot(plottingangles,normalized_firing_rate, 'LineWidth', 2);
        else
                plot(plottingangles,normalized_firing_rate, 'LineWidth', 1);
        end
        
        xlim([0 360]);
        ylim([0 0.5]);
        set(gca,'XTick',[],'YTick', [])
        
        title(cellnames(i));

        j = j + 1;
        if j>10
                j=1;
        end
end
