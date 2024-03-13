function histo_latent_overlapping(data1, data2, edgesX, edgesY, cmap1, cmap2, legendTXT1, legendTXT2, savedFileTXT, flipColor1, flipColor2, cmap_factor, output_folder)
    Nplot1 = hist3(data1','edges',{edgesX edgesY});
    Nplot1 = Nplot1';

    Nplot2 = hist3(data2','edges',{edgesX edgesY});
    Nplot2 = Nplot2';

    data1_percentage = NaN(size(Nplot1));
    overlapping_percentage = NaN(size(Nplot1));
    intersect_area = NaN(size(Nplot1));
    for i = 1:size(Nplot1, 1)
        for j = 1:size(Nplot1, 2)
            if Nplot1(i, j) == 0 & Nplot2(i, j) == 0
                continue
            end
            data1_percentage(i, j) = Nplot1(i, j) / (Nplot1(i, j) + Nplot2(i, j));
            if data1_percentage(i, j) > 0.5
                overlapping_percentage(i, j) = 1 - data1_percentage(i, j);
            else
                overlapping_percentage(i, j) = data1_percentage(i, j);
            end

            thres = 0.83;
            if data1_percentage(i, j) > thres
                intersect_area(i, j) = 1;
            elseif data1_percentage(i, j) < 1-thres
                intersect_area(i, j) = -1;
            else
                intersect_area(i, j) = 0;
            end
        end
    end
    
    figure(1)
    
    % Define the reversed gray colormap
    if flipColor1
        cmap1 = eval(['flipud(' cmap1 '(floor(length(intersect_area)*cmap_factor)+1));']);
    else
        cmap1 = eval([cmap1 '(floor(length(intersect_area)*cmap_factor));']);
        cmap1 = [[1 1 1]; cmap1];
    end
    if flipColor2
        cmap2 = eval(['flipud(' cmap2 '(floor(length(intersect_area)*cmap_factor)+1));']);
    else
        cmap2 = eval([cmap2 '(floor(length(intersect_area)*cmap_factor));']);
        cmap2 = [[1 1 1]; cmap2];
    end
    
    % Map the normalized data to RGB values using the colormap
    rgb_image = zeros([size(intersect_area) 3]);
    for i = 1:size(intersect_area, 1)
        for j = 1:size(intersect_area, 2)
            if isnan(intersect_area(i,j))
                rgb_image(i,j,:) = [1 1 1];
            elseif intersect_area(i,j) == 1
                rgb_image(i,j,:) = cmap1(floor(size(cmap1,1)/2),:);
            elseif intersect_area(i,j) == -1
                rgb_image(i,j,:) = cmap2(floor(size(cmap2,1)/2),:);
            else
                rgb_image(i,j,:) = (cmap1(floor(size(cmap1,1)/2),:) + cmap2(floor(size(cmap2,1)/2),:))/2;
            end
        end
    end
    % calculate the percentage of intersection
    intersection_percentage = sum(intersect_area == 0, 'all')/(sum(intersect_area == -1, 'all') + sum(intersect_area == 0, 'all') + sum(intersect_area == 1, 'all'));
    only1_percentage = sum(intersect_area == 1, 'all')/(sum(intersect_area == -1, 'all') + sum(intersect_area == 0, 'all') + sum(intersect_area == 1, 'all'));
    only2_percentage = sum(intersect_area == -1, 'all')/(sum(intersect_area == -1, 'all') + sum(intersect_area == 0, 'all') + sum(intersect_area == 1, 'all'));
    
    [X, Y] = meshgrid(edgesX, edgesY);
    surf(X, Y, intersect_area, rgb_image, 'EdgeColor', 'none','FaceAlpha',0.9);
    hold on
    grid off
    ylim([min(edgesY)-0.1 max(edgesY)+0.1]);
    xlim([min(edgesX)-0.1 max(edgesX)+0.1]);
    box on
    ax = gca;
    ax.LineWidth = 2;
    axis square;
    view(2)
    xlabel('z(1)','FontSize', 15);
    ylabel('z(2)','FontSize', 15);
    
    % Add legend
    % Add legend with color indications
    % h = plot(NaN, NaN, 's', 'MarkerFaceColor', cmap(floor(length(Nplot)/2), :), 'MarkerEdgeColor', cmap(floor(length(Nplot)/2), :));
    
    title(['Overlapping (' num2str(round(intersection_percentage,3)*100) '%) of the latent space'],'FontSize', 15);
    %legend(h, legendTXT, 'Location', 'Northwest','FontSize', 15);
    %legend('boxoff')
    hold off

    %%saveas(gcf, ['latent_' savedFileTXT '.jpg'], 'jpg');
    exportgraphics(gca,[output_folder filesep 'intersection_latent_' savedFileTXT '.jpg'],'Resolution',600)
    save([output_folder filesep 'intersection_latent_' savedFileTXT '.mat'], 'edgesX', 'edgesY', 'intersect_area', 'rgb_image')

    % pie chart
    % Data
    values = [intersection_percentage, only1_percentage, only2_percentage];
    labels = {[legendTXT1 '+' legendTXT2 '(' num2str(round(intersection_percentage,3)*100) '%)'], [legendTXT1 ' only' '(' num2str(round(only1_percentage, 3)*100) '%)'], [legendTXT2 ' only' '(' num2str(round(only2_percentage, 3)*100) '%)']};
    colors = {cmap2(floor(size(cmap2,1)/2),:), cmap1(floor(size(cmap1,1)/2),:), (cmap1(floor(size(cmap1,1)/2),:) + cmap2(floor(size(cmap2,1)/2),:))/2};
    
    % Pie chart
    p = pie(values, labels);
    
    % Adjust colors
    h = findobj(gca, 'Type', 'patch');
    for i = 1:numel(h)
        set(h(i), 'FaceColor', colors{i});
    end

    % Adjust line color to white
    for i = 1:2:numel(p)
        set(p(i), 'EdgeColor', 'w'); % Change 'w' to the desired line color
    end
    
    % Title
    title(['Pie chart for coverage area in the latent space']);
    
    % Legend
    %legend(labels, 'Location', 'Best');
    exportgraphics(gca,[output_folder filesep 'pie_chart' savedFileTXT '.jpg'],'Resolution',600)
end