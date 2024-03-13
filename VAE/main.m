%% MAIN
% Script for producing the figures of latent space and reconstructed data
% using pretrained variational autoencoder

clear variables %#ok<*NASGU>

folder = fileparts(which('main')); 
addpath(genpath(folder));

%% Prepare dataset
% Load dataset
data = 'fig5_CD271'; % 'fig3_CXCL12' 'fig4_CD146' 'fig5_CD271'
data_file = [folder filesep 'datasets' filesep data filesep 'lightSheetDataMLP'];
data_imported = load(data_file);
X_raw_imported = data_imported.X;
Y_raw_imported = data_imported.Y;

% Normalization
R0 = 50;
R = 0:10:200;
R = kron(R, ones(1, 5));
Factor = R0^3 ./ (R0 + R).^3;
X_raw = X_raw_imported .* Factor';
X_raw = log(X_raw + 1);

%% Load pretrained VAE
load(['pretrained_models' filesep data filesep 'pretrained_vae.mat'])

%% Load parameters setup for VAE
load(['pretrained_models' filesep data filesep 'parameters.mat'])

%% Calculate latent representations
% Training-test split
XTrain = X_raw;
XTrain(:, randomIntegers) = [];
XTest = X_raw(:, randomIntegers);
Y_raw = Y_raw_imported';
YTest = Y_raw(randomIntegers, 1);

% Calcluate latent representation for test set
dsXTest = arrayDatastore(XTest, IterationDimension=2);
dsYTest = arrayDatastore(YTest);
dsTest = combine(dsXTest, dsYTest);

mbqTest = minibatchqueue( ...
    dsTest, numOutputs, ...
    MiniBatchSize = miniBatchSize, ...
    MiniBatchFcn = @preprocessMiniBatchMLP, ...
    MiniBatchFormat= ["CB",""]);

reset(mbqTest);
[ZLatent, YLatent] = encoderPredictionsMLP(netE, mbqTest);

%% Plot latent space
outputfolder = ['replicate_' data];
mkdir(outputfolder);

idx_aged = YLatent == aged_code;
idx_young = YLatent == young_code;
aged_zlatent = ZLatent(:, idx_aged);
young_zlatent = ZLatent(:, idx_young);

histo_latent(aged_zlatent, bin_set_x, bin_set_y, 'sky', 'Aged', ['aged_' model_name], 0, cmap_factor, outputfolder)
histo_latent(young_zlatent, bin_set_x, bin_set_y, 'pink', 'Young', ['young_' model_name], 1, cmap_factor, outputfolder)
histo_latent_2surf(aged_zlatent, young_zlatent, bin_set_x, bin_set_y, 'sky', 'pink', 'Aged', 'Young', model_name, 0, 1, cmap_factor, outputfolder)
histo_latent_overlapping(aged_zlatent, young_zlatent, bin_set_x, bin_set_y, 'sky', 'pink', 'Aged', 'Young', model_name, 1, 0, 1, outputfolder)

%% Plot reconstructed data for specified regions on latent space
for i = 1:length(start_z1s)
    regions{i} = [start_z1s(i) start_z2s(i); start_z1s(i)+region_length start_z2s(i)+region_length];
    region_names{i} = ['region_' num2str(i)];
end
for i = 1:length(regions)
    cluster = regions{i};
    cluster_num = region_names{i};
    [Z_region, X_region, X_region_generated, Y_region, indexes_region] = encoderPredictionsMLP_cluster(netD,netE,mbqTest,cluster,randomIntegers);
    cellNum = length(indexes_region);
    datapoints_folder = [folder filesep 'datasets' filesep data];
    map_coord_individual(indexes_region, [cluster_num '_coordinates_' model_name], Z_region, subjects, dim, csv_base, datapoints_folder, outputfolder);
    mean_X_cluster = reconstruct_mean_input(X_region);
    mean_X_cluster_generated = reconstruct_mean_input(X_region_generated);
    plot_avg_data_distro(mean_X_cluster, mean_X_cluster_generated, 'Log(K(r) * Counts + 1)', model_name, cluster_num, cluster, celltypes, cellcolor, cellNum, outputfolder)
    close all
end