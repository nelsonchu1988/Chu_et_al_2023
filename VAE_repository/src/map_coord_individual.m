function map_coord_individual(indexes, file_name, z_coords, subs, dim, csv_base, dataset_folder, outputfolder)

    %subs = {'OP1', 'OP2','OP3', 'YP3', 'YP4', 'YP5'};
    %dim = 'CXCL12';
    z_coords = gather(z_coords);
    r = 0:10:200;
    for i = 1:length(subs)
        sub = subs{i};
        resultFile = [dataset_folder filesep sub '_dist_' dim '_' num2str(r(1)) '.mat'];
        data{i} = load(resultFile);
        counter_per_subs(i) = length(data{i}.distances);
        dataCXCL12{i} = readmatrix([dataset_folder filesep 'csv_files' filesep sub filesep sub '_' csv_base '.csv']);
    end
   
    % Find the corresponding elements for each index
    for i = 1:numel(indexes)
        index = indexes(i);
        for j = 1:numel(counter_per_subs)
            if index <= sum(counter_per_subs(1:j))
                corresponding_subs(i) = subs(j);
                coords = dataCXCL12{j};
                if j == 1
                    index_ind = index;
                else
                    index_ind = index - sum(counter_per_subs(1:j-1));
                end
                coord_x(i) = coords(index_ind, 1);
                coord_y(i) = coords(index_ind, 2);
                coord_z(i) = coords(index_ind, 3);
                break;
            end
        end
    end
    % % % double check the latent space
    % % % Nplot = hist3(z_coords','edges',{-3:.1:3 -3:.1:3});
    % % % Nplot = Nplot';
    % % % [X, Y] = meshgrid(-3:.1:3, -3:.1:3);
    % % % surf(X, Y, Nplot, 'EdgeColor', 'none');
    % % % view(2)
    latent_z1 = z_coords(1, :);
    latent_z2 = z_coords(2, :);
    T = table(corresponding_subs', coord_x', coord_y', coord_z', latent_z1', latent_z2', 'VariableNames', ["Individual", "X", "Y", "Z", "latent_Z1", "latent_Z2"]);
    writetable(T, [outputfolder filesep file_name '.xlsx']);
end