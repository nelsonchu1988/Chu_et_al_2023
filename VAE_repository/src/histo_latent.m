function histo_latent(data, edgesX, edgesY, cmap, legendTXT, savedFileTXT, flipColor, cmap_factor, output_folder)
    Nplot = hist3(data','edges',{edgesX edgesY});
    Nplot = Nplot';
    figure(1)
    
    % Define the reversed gray colormap
    if flipColor
        cmap = eval(['flipud(' cmap '(floor(length(Nplot)*cmap_factor)+1));']);
    else
        cmap = eval([cmap '(floor(length(Nplot)*cmap_factor));']);
        cmap = [[1 1 1]; cmap];
    end
    
    % Map the normalized data to RGB values using the colormap
    rgb_image = ind2rgb(Nplot, cmap);
    
    [X, Y] = meshgrid(edgesX, edgesY);
    surf(X, Y, Nplot, rgb_image, 'EdgeColor', 'none');
    hold on
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
%    h = plot(NaN, NaN, 's', 'MarkerFaceColor', cmap(floor(length(Nplot)/2), :), 'MarkerEdgeColor', cmap(floor(length(Nplot)/2), :));
    
    title('Visualization of the latent space','FontSize', 15);
    %legend(h, legendTXT, 'Location', 'Northwest','FontSize', 15);
    %legend('boxoff')
    hold off

    %%saveas(gcf, ['latent_' savedFileTXT '.jpg'], 'jpg');
    exportgraphics(gca,[output_folder filesep 'latent_' savedFileTXT '.jpg'],'Resolution',600)
    save([output_folder filesep 'latent_' savedFileTXT '.mat'], 'X', 'Y', 'Nplot', 'rgb_image')
end