function plot_avg_data_distro(mean_X_cluster, mean_X_cluster_generated, y_label, model_name, cluster_num, cluster_cood, output_folder)
figure; hold on; 
x = 0:10:200;
% % % scatter(x, (mean_X_cluster_generated(1, :)), 'xred', 'LineWidth', 2, 'MarkerEdgeAlpha', 0.2);
% % % scatter(x, (mean_X_cluster_generated(2, :)), '+blue', 'LineWidth', 2, 'MarkerEdgeAlpha', 0.2);
% % % scatter(x, (mean_X_cluster_generated(3, :)), '<yellow', 'LineWidth', 2, 'MarkerEdgeAlpha', 0.2);
% % % scatter(x, (mean_X_cluster_generated(4, :)), '>green', 'LineWidth', 2, 'MarkerEdgeAlpha', 0.2);
% % % scatter(x, (mean_X_cluster_generated(5, :)), '*', 'color', [107 76 154]./255, 'LineWidth', 2, 'MarkerEdgeAlpha', 0.2);

plot(x, (mean_X_cluster_generated(1, :)), 'red', 'LineWidth', 2, 'LineStyle', '--');
plot(x, (mean_X_cluster_generated(2, :)), 'blue', 'LineWidth', 2, 'LineStyle', '--');
plot(x, (mean_X_cluster_generated(3, :)), 'yellow', 'LineWidth', 2, 'LineStyle', '--');
plot(x, (mean_X_cluster_generated(4, :)), 'green', 'LineWidth', 2, 'LineStyle', '--');
plot(x, (mean_X_cluster_generated(5, :)), 'color', [107 76 154]./255, 'LineWidth', 2, 'LineStyle', '--');

% % % plot(x, (mean_X_cluster(1, :)), 'red', 'LineWidth', 2);
% % % plot(x, (mean_X_cluster(2, :)), 'blue', 'LineWidth', 2);
% % % plot(x, (mean_X_cluster(3, :)), 'yellow', 'LineWidth', 2);
% % % plot(x, (mean_X_cluster(4, :)), 'green', 'LineWidth', 2);
% % % plot(x, (mean_X_cluster(5, :)), 'color', [107 76 154]./255, 'LineWidth', 2);


%set(gca, 'YScale', 'log')
xticks(0:40:200);
ylim([0 2])

yticks(0:0.5:2);
% Customize the plot
ax = gca;
ax.LineWidth = 2;
set(gca,'fontsize',15)

% Customize the plot
box on;
xlabel('Radius (\mum)', 'FontSize', 15); % Include the unit "μm" in the label
ylabel(y_label, 'FontSize', 15);

% Set axis limits to have an equal aspect ratio
axis square;

% Add grid
%grid on;

% Title and legend
title(['region z(1): [' num2str(cluster_cood(1)) ', ' num2str(cluster_cood(2)) '], z(2): [' num2str(cluster_cood(3)) ', ' num2str(cluster_cood(4)) ']'],'FontSize', 15);
%legend('CXCL Generated', 'BONE Generated', 'CD31 Generated', 'YOPRO Generated', 'PER Generated', 'Location', 'northwest','FontSize', 6);
%legend('boxoff')

hold off

exportgraphics(gca, [output_folder filesep 'avg_data_distro_' cluster_num '_' model_name '.jpg'],'Resolution',600);
save([output_folder filesep 'avg_data_distro_' cluster_num '_' model_name '.mat'], 'mean_X_cluster', 'mean_X_cluster_generated')
end