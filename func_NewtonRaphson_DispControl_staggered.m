function [J_uu, DRdu, dofs, Res_uu_F_norm, Global_Residual_norm, iteration, history_var_mat, gausspoints_prop_mat, nodes_prop_mat, Reactions_y, Res_ee_norm, time_counter_val] = func_NewtonRaphson_DispControl_staggered(dofs,fixnodes_applied,increment,Delastic,history_var_mat,num_elem_at_node,n_hood,weights,IsProj,model_name,loadfactor,loadfactor_lim,Modified_NR_flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ==================== NEWTON RAPHSON ANALYSIS METHOD =====================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize the time counter 
time_counter_val = 0;

% Include global variables
func_include_flags;

if IsProj == 0
    disp(strcat("ndof = ", num2str(ndof)))
end

% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOLVER ROUTINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
if IsProj == 0

    % Create empty entries for element/nodal properties 
    gausspoints_prop_mat = [];
    nodes_prop_mat       = [];   
    DRdu                 = [];
    Res_uu_F_norm        = [];
    Res_ee_norm          = [];
    
    % Create copies of the kappa matrix from the end of the previous increment
    history_var_mat_previousinc = history_var_mat;
    
    % Initiate the iteration counter within each loading step
    iteration = 1;
    
    % Evaluate stiffness matrix [K] and history variable matrix at the start of the increment
    [J_uu, ~, ~, history_var_mat, ~, ~, ~, ~, time_counter_val] = func_globalstiffness_staggered_u(dofs,Delastic,model_name,increment,iteration,history_var_mat_previousinc,num_elem_at_node,n_hood,weights,IsProj,loadfactor,loadfactor_lim,time_counter_val,Modified_NR_flag);

    % Partition K and DRdu
    if TangentID == 1
        [~, ~, ~, J_uu_F, ~, ~, ~, ~] = func_partitionK_FOR_REACTIONS(fixnodes_applied, J_uu); 
    else
        disp("Check your flags! - func_NewtonRaphson_DispControl_staggered")
    end

    % Separate displacements from strains in the dof vector:
    dofs_uu_x = dofs(1:3:length(dofs)-2);
    dofs_uu_y = dofs(2:3:length(dofs)-1);            
    dofs_uu   = reshape([dofs_uu_x dofs_uu_y]',[],1);
    dofs_ee   = dofs(3:3:length(dofs));
    [~, dofs_E, dofs_F, ~, ~, listofnodes_ebc_global, listofnodes_nbc_global] = func_partitiond_fixnodes(fixnodes_applied,dofs);
    [~, dofs_uu_E, dofs_uu_F, ~, ~, listofnodes_ebc_uu, listofnodes_nbc_uu]   = func_partitiond_fixnodes_FOR_REACTIONS(fixnodes_applied,dofs_uu);

    % Assemble the global dofs
    dofs(listofnodes_ebc_global) = dofs_E; 
    dofs(listofnodes_nbc_global) = dofs_F;
    
    % Calculate and partition the residual force vector (based on internal stresses) after delta_dofs_E has been applied
    % disp("Next call of globalstiffness: Obtain Residual (after just dE has been applied)")
    [~, ~, Res_uu, ~, ~, ~, ~, ~, time_counter_val] = func_globalstiffness_staggered_u(dofs,Delastic,model_name,increment,iteration,history_var_mat_previousinc,num_elem_at_node,n_hood,weights,IsProj,loadfactor,loadfactor_lim,time_counter_val,Modified_NR_flag);
    [~, ~, Res_uu_F] = func_partitiond_fixnodes_FOR_REACTIONS(fixnodes_applied,Res_uu); 
        
    % The while loop runs until either the error is smaller than a specific
    % threshold (convergence is achieved) or the number of iterations within 
    % a loading step is exceeded (computational cost is high)
    while true
    
        % Calculate the displacements increment at the natural boundary
        if TangentID == 1
            delta_dofs_uu_F = - J_uu_F \ Res_uu_F;     % Analytical calculation (using K)
        else
            disp("Check your flags! - func_NewtonRaphson_DispControl_staggered")
        end
        
        % Update the displacements at the natural boundary 
        dofs_uu_F = dofs_uu_F + delta_dofs_uu_F;
        
        % Assemble the displacement vector
        dofs_uu(listofnodes_ebc_uu) = dofs_uu_E;
        dofs_uu(listofnodes_nbc_uu) = dofs_uu_F;

        % Assemble the global dofs vector
        dofs(1:3:length(dofs)-2) = dofs_uu(1:2:length(dofs_uu)-1);
        dofs(2:3:length(dofs)-1) = dofs_uu(2:2:length(dofs_uu));
        dofs(3:3:length(dofs))   = dofs_ee;

        % Compute tangent and residual based on the nonlocal strain
        [J_ee, ~, Res_ee, ~, ~, ~, ~, ~, time_counter_val] = func_globalstiffness_staggered_e(dofs,Delastic,model_name,increment,iteration,history_var_mat_previousinc,num_elem_at_node,n_hood,weights,IsProj,loadfactor,loadfactor_lim,time_counter_val,Modified_NR_flag);

        % Calculate the nonlocal strain increment
        if TangentID == 1
            delta_dofs_ee = - J_ee \ Res_ee;     % Analytical calculation (using K)
        else
            disp("Check your flags! - func_NewtonRaphson_DispControl_staggered")
        end

        % Update the nonlocal strain dofs
        dofs_ee = dofs_ee + delta_dofs_ee;
        dofs(3:3:length(dofs)) = dofs_ee;

        % Evaluate the new stiffness matrix, residual force vector (based on internal stresses) and history variable matrix
        %disp("Next call of globalstiffness: Obtain K, history variable and Res (inside the while loop)")
        [J_uu, ~, Res_uu, history_var_mat, ~, ~, ~, ~, time_counter_val] = func_globalstiffness_staggered_u(dofs,Delastic,model_name,increment,iteration,history_var_mat_previousinc,num_elem_at_node,n_hood,weights,IsProj,loadfactor,loadfactor_lim,time_counter_val,Modified_NR_flag);

        % Partition the stiffness matrix [K]/[DRdu] and residual force vector {Res_F}
        [~, ~, ~, J_uu_F, ~, ~, ~, ~]  = func_partitionK_FOR_REACTIONS(fixnodes_applied,J_uu);
        [~, ~, Res_uu_F]               = func_partitiond_fixnodes_FOR_REACTIONS(fixnodes_applied,Res_uu); 
        
        % Calculate the norm of the global residual
        Global_Residual_norm(iteration) = norm([delta_dofs_uu_F; delta_dofs_ee],2);
       
        % Keep track of the residual at the first iteration
        if iteration == 1 
            Global_Residual_first = Global_Residual_norm(iteration);
        end
        
        % If the residual error is less than the convergence threshold then the
        % solution is converged, otherwise proceed to the next iteration 
        if (Global_Residual_norm(iteration) / Global_Residual_first < 10^-4) || (iteration == max_accept_iter)
            
            % -------------------------------------------------------------
            % Compute the elastic stiffness matrix K = (1-D) * B' * C * B
            [~, ~, ~, ~, ~, ~, ~, K, time_counter_val] = func_globalstiffness(dofs,Delastic,model_name,increment,iteration,history_var_mat_previousinc,num_elem_at_node,n_hood,weights,IsProj,loadfactor,loadfactor_lim,time_counter_val,Modified_NR_flag);           
            
            % Partition the elastic stiffness matrix based on the essential
            % and natural boundary
            [K_E, K_EF, ~, ~, ~, ~, ~, ~]   = func_partitionK_FOR_REACTIONS(fixnodes_applied,K);

            if ndof == 3
                U_x = dofs(1:3:length(dofs)-2);
                U_y = dofs(2:3:length(dofs)-1);            
                U   = reshape([U_x U_y]',[],1);
            elseif ndof == 2
                U = dofs;
            end
            
            % Partition the degrees of freedom vector based on the essential 
            % and natural boundary
            [~, dofs_U_E, dofs_U_F, ~, ~, ~, ~] = func_partitiond_fixnodes_FOR_REACTIONS(fixnodes_applied, U);
            
            % Compute the reaction forces
            f_ext_E = K_E * dofs_U_E + K_EF * dofs_U_F; 
            Reactions_y = sum(abs(f_ext_E(fixnodes_applied(2,:) == 2 & fixnodes_applied(3,:) == 0)))

            
            Sum_Fx = sum(f_ext_E(fixnodes_applied(2,:) == 1))
            Sum_Fy = sum(f_ext_E(fixnodes_applied(2,:) == 2))

            break

        else
            iteration = iteration + 1;
        end 
        
    end


elseif IsProj == 1
    % ---------------------------------------------------------------------
    %%%%%%%%%%%%%%%%%%%%% PLOTTING PROPERTIES ROUTINE %%%%%%%%%%%%%%%%%%%%%
    % ---------------------------------------------------------------------
    
    % Create empty entries for solver variables
    J_uu                    = [];
    DRdu                    = [];
    Res_uu_F_norm           = []; 
    Global_Residual_norm    = []; 
    iteration               = [];  
    Reactions_y             = [];
    Res_ee_norm             = [];
    
    % Create copies of the kappa matrix from the end of the previous increment
    history_var_mat_previousinc = history_var_mat;
    
    % Compute element/nodal properties
    [~, ~, ~, ~, gausspoints_prop_mat, nodes_prop_mat, Res_ee_norm, ~] = func_globalstiffness(dofs,Delastic,model_name,increment,iteration,history_var_mat_previousinc,num_elem_at_node,n_hood,weights,IsProj,loadfactor,loadfactor_lim,time_counter_val,Modified_NR_flag);    


else 

    disp("Check your IsProj variable")

end




end


