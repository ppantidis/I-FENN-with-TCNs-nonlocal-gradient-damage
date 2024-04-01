clc
clear all
close all
format compact

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TrainingData_StrA100_GP_features_lf_time        = zeros(3600,5,160);
TrainingData_StrA100_Nodeslrb_features_lf_time  = zeros(82,5,160);
TrainingData_StrA100_Nodesbtb_features_lf_time  = zeros(82,5,160);
Laplacian_StrA100_GP                            = zeros(3600,160);

% Load all .mat files
for i = 1:160
    
    % Load file
    load(strcat('Data_FEM_StrA100_Quad_Nonlocal_Gradient_inc_',num2str(i),'.mat'));

    % ---------------------------------------------------------------------
    % Data at Gauss Points
    TrainingData_StrA100_GP_features_lf_time(:,:,i) = data_GP;

    % Data at Left-right boundary nodes
    TrainingData_StrA100_Nodeslrb_features_lf_time(:,:,i) = data_nodes_lrb_StrA100_Quad;

    % Data at Bottom-Top boundary nodes
    TrainingData_StrA100_Nodesbtb_features_lf_time(:,:,i) = data_nodes_btb_StrA100_Quad;

    % ---------------------------------------------------------------------
    % Laplacian terms at Gauss Points
    Partial_nonlocal_physical_transp = Partial_nonlocal_physical';
    Laplacian_StrA100_GP(:,i) = Partial_nonlocal_physical_transp(:,1) + Partial_nonlocal_physical_transp(:,2);
    

end

% Re-arrange the tensors
TrainingData_StrA100_GP_time_features_lf       = permute(TrainingData_StrA100_GP_features_lf_time,[1 3 2]);
TrainingData_StrA100_Nodeslrb_time_features_lf = permute(TrainingData_StrA100_Nodeslrb_features_lf_time,[1 3 2]);
TrainingData_StrA100_Nodesbtb_time_features_lf = permute(TrainingData_StrA100_Nodesbtb_features_lf_time,[1 3 2]);

save TrainingData_StrA100_GP_time_features_lf TrainingData_StrA100_GP_time_features_lf
save TrainingData_StrA100_Nodeslrb_time_features_lf TrainingData_StrA100_Nodeslrb_time_features_lf
save TrainingData_StrA100_Nodesbtb_time_features_lf TrainingData_StrA100_Nodesbtb_time_features_lf
save Laplacian_StrA100_GP Laplacian_StrA100_GP

whos








