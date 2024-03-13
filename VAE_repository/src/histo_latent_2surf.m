function histo_latent_2surf(data1, data2, edgesX, edgesY, cmap1, cmap2, legendTXT1, legendTXT2, savedFileTXT, flipColor1, flipColor2, cmap_factor, output_folder)
    Nplot = hist3(data1','edges',{edgesX edgesY});
    Nplot = Nplot';
    
    figure(1)
    
    % Define the reversed gray colormap
    if flipColor1
        cmap1 = eval(['flipud(' cmap1 '(floor(length(Nplot)*cmap_factor)+1));']);
    else
        cmap1 = eval([cmap1 '(floor(length(Nplot)*cmap_factor));']);
        cmap1 = [[1 1 1]; cmap1];
    end
    
    % Map the normalized data to RGB values using the colormap
    rgb_image = ind2rgb(Nplot, cmap1);
    
    [X, Y] = meshgrid(edgesX, edgesY);
    surf(X, Y, Nplot, rgb_image, 'EdgeColor', 'none');
    hold on
    
    Nplot = hist3(data2','edges',{edgesX edgesY});
    Nplot = Nplot';
    % Define the reversed gray colormap
    if flipColor2
        cmap2 = eval(['flipud(' cmap2 '(floor(length(Nplot)*cmap_factor)+1));']);
    else
        cmap2 = eval([cmap2 '(floor(length(Nplot)*cmap_factor));']);
        cmap2 = [[1 1 1]; cmap2];
    end
    
    % Map the normalized data to RGB values using the colormap
    rgb_image = ind2rgb(Nplot, cmap2);
    
    [X, Y] = meshgrid(edgesX, edgesY);
    surf(X, Y, Nplot, rgb_image, 'EdgeColor', 'none');
    

    grid off
    ylim([min(edgesY)-0.1 max(edgesY)+0.1]);
    xlim([min(edgesX)-0.1 max(edgesX)+0.1]);
    box on
    axis square;
    view(2)
    xlabel('z(1)','FontSize', 15);
    ylabel('z(2)','FontSize', 15);
    ax = gca;
    ax.LineWidth = 2;
    % Add legend
    % Add legend with color indications
%    h = zeros(2, 1);
%    h(1) = plot(NaN, NaN, 's', 'MarkerFaceColor', cmap1(floor(length(Nplot)/2), :), 'MarkerEdgeColor', cmap1(floor(length(Nplot)/2), :));
%    h(2) = plot(NaN, NaN, 's', 'MarkerFaceColor', cmap2(floor(length(Nplot)/2), :), 'MarkerEdgeColor', cmap2(floor(length(Nplot)/2), :));
    
    title('Visualization of the latent space','FontSize', 15);
    %legend(h, legendTXT1, legendTXT2, 'Location', 'Northwest','FontSize', 15);
    %legend('boxoff')
    hold off

    % % saveas(gcf, ['latent_' savedFileTXT '.jpg'], 'jpg');
    exportgraphics(gca,[output_folder filesep 'latent_' savedFileTXT '.jpg'],'Resolution',600)
    save([output_folder filesep 'latent_' savedFileTXT '.mat'], 'X', 'Y', 'Nplot', 'rgb_image')
end