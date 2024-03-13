

function reconstructed_mean_x = reconstruct_mean_input(x_points)
reconstructed_mean_x = [];
% CXCL -> idx 1
idx_cxcl = 1:5:105;
reconstructed_mean_x(1, :) = mean(x_points(idx_cxcl, :), 2);

% Bone -> idx 2
idx_bone = 2:5:105;
reconstructed_mean_x(2, :) = mean(x_points(idx_bone, :), 2);

% Cd31 -> idx 3
idx_cd31 = 3:5:105;
reconstructed_mean_x(3, :) = mean(x_points(idx_cd31, :), 2);

% yopro -> idx 4
idx_yopro = 4:5:105;
reconstructed_mean_x(4, :) = mean(x_points(idx_yopro, :), 2);

% peri -> idx 5
idx_peri = 5:5:105;
reconstructed_mean_x(5, :) = mean(x_points(idx_peri, :), 2);
end