function [] = func_savingresults(xcoord_GP_mat,ycoord_GP_mat,g,loadfactor,gausspoints_prop_mat,increment,model_name,solver_string,tangent_string,nnodes,lelem,nodes_prop_mat,Res_F_F_norm,Res_u_norm,Data_IFENN_GP,Res_e_norm,Reactions_mat,coords,Partial_nonlocal_physical,time_counter,Modified_NR_flag)
% ======================= SAVING RESULTS ==================================

if Modified_NR_flag == 0
    flag_NR = "FNR";
elseif Modified_NR_flag == 1
    flag_NR = "MNR";
else
    disp('Check your Modified NR flag')
end

% -------------------------------------------------------------------------
% Gauss Points
xcoord_GP_vec       = reshape(xcoord_GP_mat',1,[ ])';
ycoord_GP_vec       = reshape(ycoord_GP_mat',1,[ ])';
loadfactor_GP_vec   = ones(length(xcoord_GP_vec),1) * loadfactor;
data_GP             = [xcoord_GP_vec ycoord_GP_vec gausspoints_prop_mat(:,8:9) loadfactor_GP_vec gausspoints_prop_mat(:,1)];

if model_name == "StrA100_Quad"

    % -------------------------------------------------------------------------
    % Left/right boundary 
    xcoord_nodes_lrb_vec        = [zeros(length(0:2.5:100),1); zeros(length(0:2.5:100),1)];
    ycoord_nodes_lrb_vec        = [0:2.5:100 0:2.5:100]';
    ID_nodes_lrb_vec            = [find(coords(1,:)==0) find(coords(1,:)==100)]';
    nodes_lrb_prop_mat          = nodes_prop_mat(ID_nodes_lrb_vec,:);
    loadfactor_nodes_lrb_vec    = ones(length(xcoord_nodes_lrb_vec),1) * loadfactor;
    data_nodes_lrb_StrA100_Quad = [xcoord_nodes_lrb_vec ycoord_nodes_lrb_vec nodes_lrb_prop_mat(:,8:9) loadfactor_nodes_lrb_vec];
    
    
    % -------------------------------------------------------------------------
    % Top/bottom boundary
    xcoord_nodes_btb_vec        = [0:2.5:100 0:2.5:100]';
    ycoord_nodes_btb_vec        = [zeros(length(0:2.5:100),1); zeros(length(0:2.5:100),1)];
    ID_nodes_btb_vec            = [find(coords(2,:)==0) find(coords(2,:)==100)]';
    nodes_btb_prop_mat          = nodes_prop_mat(ID_nodes_btb_vec,:);
    loadfactor_nodes_btb_vec    = ones(length(xcoord_nodes_btb_vec),1) * loadfactor;
    data_nodes_btb_StrA100_Quad = [xcoord_nodes_btb_vec ycoord_nodes_btb_vec nodes_btb_prop_mat(:,8:9) loadfactor_nodes_btb_vec];

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVING IN .MAT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Save the results for this increment
matfilename = strcat('Data_', model_name, '_', solver_string, '_', flag_NR, '_inc_', num2str(increment), '.mat');

if model_name == "StrA100_Quad"
    % Save data for Gauss Points and boundary nodes
    save(matfilename, 'data_GP', 'data_nodes_lrb_StrA100_Quad', 'data_nodes_btb_StrA100_Quad', 'Res_F_F_norm', 'Res_u_norm', 'Res_e_norm', 'Reactions_mat', 'Partial_nonlocal_physical', 'time_counter'); 
else
    % Save data for Gauss Points
    save(matfilename, 'data_GP', 'Res_F_F_norm', 'Res_u_norm', 'Res_e_norm', 'Reactions_mat', 'Partial_nonlocal_physical', 'time_counter'); 
end

% % -------------------------------------------------------------------------
% % Update the "Data_IFENN_GP" tensor which is used as input for the next TCN
% % prediction
% Data_IFENN_GP(:,increment,:) = [xcoord_GP_vec ycoord_GP_vec gausspoints_prop_mat(:,8) gausspoints_prop_mat(:,9) loadfactor_GP_vec];
% save Data_IFENN_GP Data_IFENN_GP


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVING IN EXCEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% cheader = ["xcoord" "ycoord" "g" "damage" "sigma_xx" "sigma_yy" "tau_xy" "e_xx" "e_yy" "e_xy" "elocal" "enonlocal_pred"];
% commaHeader = [cheader;repmat({','},1,numel(cheader))]; % insert commas
% commaHeader = commaHeader(:)';
% textHeader  = cell2mat(commaHeader);                    % cHeader in text with commas
% 
% % -------------------------------------------------------------------------
% % Save data for GP
% fid1 = fopen(strcat("data_GP_" + model_name + "_" + solver_string + "_" + tangent_string + "_inc_" + int2str(increment) + ".csv"),'w'); 
% fprintf(fid1,'%s\n',textHeader);
% fclose(fid1);
% dlmwrite(strcat("data_GP_" + model_name + "_" + solver_string + "_" + tangent_string + "_inc_" + int2str(increment) + ".csv"), data_GP_StrA100,'-append', 'precision', 16);
% 
% % -------------------------------------------------------------------------
% % Save data for left-right boundary nodes
% fid2 = fopen(strcat("data_nodes_lrb_" + model_name + "_" + solver_string + "_" + tangent_string + "_inc_" + int2str(increment) + ".csv"),'w'); 
% fprintf(fid2,'%s\n',textHeader);
% fclose(fid2);
% dlmwrite(strcat("data_nodes_lrb_" + model_name + "_" + solver_string + "_" + tangent_string + "_inc_" + int2str(increment) + ".csv"), data_nodes_lrb_StrA100,'-append', 'precision', 16);
% 
% % -------------------------------------------------------------------------
% % Save data for top-bottom boundary nodes
% fid3 = fopen(strcat("data_nodes_btb_" + model_name + "_" + solver_string + "_" + tangent_string + "_inc_" + int2str(increment) + ".csv"),'w'); 
% fprintf(fid3,'%s\n',textHeader);
% fclose(fid3);
% dlmwrite(strcat("data_nodes_btb_" + model_name + "_" + solver_string + "_" + tangent_string + "_inc_" + int2str(increment) + ".csv"), data_nodes_btb_StrA100,'-append', 'precision', 16);


end









