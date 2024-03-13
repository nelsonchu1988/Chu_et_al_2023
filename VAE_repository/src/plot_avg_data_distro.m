function plot_avg_data_distro(mean_X_cluster, mean_X_cluster_generated, y_label, model_name, cluster_num, cluster_cood, celltypes, cellcolor, cellNum, outputfolder)
figure; hold on; 
x = 0:10:200;
% % % scatter(x, (mean_X_cluster_generated(1, :)), ['x' cellcolor{1}], 'LineWidth', 2, 'MarkerEdgeAlpha', 0.2);
% % % scatter(x, (mean_X_cluster_generated(2, :)), ['+' cellcolor{2}], 'LineWidth', 2, 'MarkerEdgeAlpha', 0.2);
% % % scatter(x, (mean_X_cluster_generated(3, :)), ['<' cellcolor{3}], 'LineWidth', 2, 'MarkerEdgeAlpha', 0.2);
% % % scatter(x, (mean_X_cluster_generated(4, :)), ['>' cellcolor{4}], 'LineWidth', 2, 'MarkerEdgeAlpha', 0.2);
% % % scatter(x, (mean_X_cluster_generated(5, :)), ['*' cellcolor{5}], 'LineWidth', 2, 'MarkerEdgeAlpha', 0.2);
% % % 
% % % plot(x, (mean_X_cluster(1, :)), cellcolor{1}, 'LineWidth', 2);
% % % plot(x, (mean_X_cluster(2, :)), cellcolor{2}, 'LineWidth', 2);
% % % plot(x, (mean_X_cluster(3, :)), cellcolor{3}, 'LineWidth', 2);
% % % plot(x, (mean_X_cluster(4, :)), cellcolor{4}, 'LineWidth', 2);
% % % plot(x, (mean_X_cluster(5, :)), cellcolor{5}, 'LineWidth', 2);

plot(x, (mean_X_cluster_generated(1, :)), cellcolor{1}, 'LineWidth', 2, 'LineStyle', '--');
plot(x, (mean_X_cluster_generated(2, :)), cellcolor{2}, 'LineWidth', 2, 'LineStyle', '--');
plot(x, (mean_X_cluster_generated(3, :)), cellcolor{3}, 'LineWidth', 2, 'LineStyle', '--');
plot(x, (mean_X_cluster_generated(4, :)), cellcolor{4}, 'LineWidth', 2, 'LineStyle', '--');
if isnumeric(cellcolor{5})
    plot(x, (mean_X_cluster_generated(5, :)), 'color', cellcolor{5}, 'LineWidth', 2, 'LineStyle', '--');
else
    plot(x, (mean_X_cluster_generated(5, :)), cellcolor{5}, 'LineWidth', 2, 'LineStyle', '--');
end
%set(gca, 'YScale', 'log')
xticks(0:40:200);
yticks(0:0.5:2);
ylim([0 2])

% Customize the plot
box on;
xlabel('Radius (\mum)', 'FontSize', 15); % Include the unit "μm" in the label
ylabel(y_label, 'FontSize', 15);

% Set axis limits to have an equal aspect ratio
axis square;

% Add grid
grid off;

% Title and legend
title(['cluster z(1): [' num2str(cluster_cood(1)) ', ' num2str(cluster_cood(2)) '], z(2): [' num2str(cluster_cood(3)) ', ' num2str(cluster_cood(4)) ']'],'FontSize', 15);
%legend(celltypes{1}, celltypes{2}, celltypes{3}, celltypes{4}, celltypes{5}, 'Location', 'northwest','FontSize', 6);
%legend('boxoff')

hold off
exportgraphics(gca, [outputfolder filesep 'avg_data_distro_' cluster_num '_' model_name '.jpg'],'Resolution',600);
save([outputfolder filesep 'avg_data_distro_' cluster_num '_' model_name '.mat'], 'mean_X_cluster', 'mean_X_cluster_generated')

end