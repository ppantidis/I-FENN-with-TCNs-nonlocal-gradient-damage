function [J, DRdu, dofs, Res_F_F_norm, Res_u_norm, iteration, history_var_mat, gausspoints_prop_mat, nodes_prop_mat, Reactions_y, Res_e_norm, time_counter_val, IFENN_flag] = func_NewtonRaphson_DispControl(dofs,fixnodes_applied,increment,Delastic,history_var_mat,num_elem_at_node,n_hood,weights,IsProj,model_name,loadfactor,loadfactor_lim,Modified_NR_flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ==================== NEWTON RAPHSON ANALYSIS METHOD =====================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize the time counter 
time_counter_val = 0;

% Include global variables
func_include_flags;

% Switch the solver
if loadfactor > loadfactor_lim
    ndof = 2;
    IFENN_flag = 1;
else
    IFENN_flag = 0;
end

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

    % Create copies of the kappa matrix from the end of the previous increment
    history_var_mat_previousinc = history_var_mat;
    
    % Initiate the iteration counter within each loading step
    iteration = 1;
    
    % Evaluate stiffness matrix [K] and history variable matrix at the start of the increment
    [J, DRdu, ~, history_var_mat, ~, ~, ~, ~, time_counter_val] = func_globalstiffness(dofs,Delastic,model_name,increment,iteration,history_var_mat_previousinc,num_elem_at_node,n_hood,weights,IsProj,loadfactor,loadfactor_lim,time_counter_val,Modified_NR_flag);

    % Partition K and DRdu
    if TangentID == 1
        [~, ~, ~, J_F, ~, ~, listofnodes_ebc, listofnodes_nbc] = func_partitionK(fixnodes_applied, J); 
    elseif TangentID == 2
        [~, ~, ~, DRdu_F, ~, ~, listofnodes_ebc, listofnodes_nbc] = func_partitionK(fixnodes_applied, DRdu);
    else
        disp("Check your choice of solver!")
    end
    
    % Calculate the fixed displacement increment and partition dofs
    [~, dofs_E, dofs_F] = func_partitiond_fixnodes(fixnodes_applied,dofs);
    
    % Assemble the global dofs
    dofs(listofnodes_ebc) = dofs_E; 
    dofs(listofnodes_nbc) = dofs_F;
    
    % Calculate and partition the residual force vector (based on internal stresses) after delta_dofs_E has been applied
    % disp("Next call of globalstiffness: Obtain Residual (after just dE has been applied)")
    [~, ~, Res_F, ~, ~, ~, Res_e_first, ~, time_counter_val] = func_globalstiffness(dofs,Delastic,model_name,increment,iteration,history_var_mat_previousinc,num_elem_at_node,n_hood,weights,IsProj,loadfactor,loadfactor_lim,time_counter_val,Modified_NR_flag);
    [~, ~, Res_F_F] = func_partitionf(fixnodes_applied,Res_F); 
    
    % Compute and store the norm of the residual at the free boundary 
    Res_F_F_norm_first      = norm(Res_F_F,2);
    Res_F_F_norm(iteration) = Res_F_F_norm_first;
    
    % ---------------------------------------------------------------------
    Res_e_norm_first      = norm(Res_e_first,2);
    Res_e_norm(iteration) = Res_e_norm_first;
    % ---------------------------------------------------------------------

    % The while loop runs until either the error is smaller than a specific
    % threshold (convergence is achieved) or the number of iterations within 
    % a loading step is exceeded (computational cost is high)
    while true
                
        % Calculate the displacements increment at the natural boundary
        if TangentID == 1
            delta_dofs_F = - J_F \ Res_F_F;     % Analytical calculation (using K)
        elseif TangentID == 2
            delta_dofs_F = - DRdu_F \ Res_F_F;  % Numerical approximation (using DRdu) 
        else
            disp("Check your choice of solver!")
        end
        
        % Update the displacements at the natural boundary 
        dofs_F = dofs_F + delta_dofs_F;
    
        % Assemble the global dofs vector
        dofs(listofnodes_ebc) = dofs_E;
        dofs(listofnodes_nbc) = dofs_F;
       
        % Evaluate the new stiffness matrix, residual force vector (based on internal stresses) and history variable matrix
        %disp("Next call of globalstiffness: Obtain K, history variable and Res (inside the while loop)")
        [J, DRdu, Res_F, history_var_mat, ~, ~, Res_e, ~, time_counter_val] = func_globalstiffness(dofs,Delastic,model_name,increment,iteration,history_var_mat_previousinc,num_elem_at_node,n_hood,weights,IsProj,loadfactor,loadfactor_lim,time_counter_val,Modified_NR_flag);

        % Partition the stiffness matrix [K]/[DRdu] and residual force vector {Res_F}
        [~, ~, ~, J_F, ~, ~, ~, ~]     = func_partitionK(fixnodes_applied,J);
        [~, ~, ~, DRdu_F, ~, ~, ~, ~]  = func_partitionK(fixnodes_applied,DRdu);
        [~, ~, Res_F_F]                = func_partitionf(fixnodes_applied,Res_F); 
        
        % Calculate the norm of the residuals based on: 
        % a) internal stresses (Res_F_F_norm) 
        % b) displacements (delta_dofs_F)
        % c) non-local strains (Res_e_norm)
        Res_F_F_norm(iteration+1)   = norm(Res_F_F,2);
        Res_u_norm(iteration)       = norm(delta_dofs_F,2); 
        Res_e_norm(iteration+1)     = norm(Res_e,2);
        
        % Keep track of the residual at the first iteration
        if iteration == 1 
            Res_u_norm_first = Res_u_norm(iteration);
        end
        
        % If the residual error is less than the convergence threshold then the
        % solution is converged, otherwise proceed to the next iteration 
        if (Res_u_norm(iteration) / Res_u_norm_first < 10^-6) || (iteration == max_accept_iter)
            
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
%             Reactions_y = sum(abs(f_ext_E(fixnodes_applied(2,:) == 2 & fixnodes_applied(3,:) == 0)))
            Reactions_y = sum(abs(f_ext_E(fixnodes_applied(2,:) == 1 & fixnodes_applied(3,:) ~= 0)))
           
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
    J            = [];
    DRdu         = [];
    Res_F_F_norm = []; 
    Res_u_norm   = []; 
    iteration    = [];  
    Reactions_y  = [];
    
    % Create copies of the kappa matrix from the end of the previous increment
    history_var_mat_previousinc = history_var_mat;
    
    % Compute element/nodal properties
    [~, ~, ~, ~, gausspoints_prop_mat, nodes_prop_mat, Res_e_norm, ~] = func_globalstiffness(dofs,Delastic,model_name,increment,iteration,history_var_mat_previousinc,num_elem_at_node,n_hood,weights,IsProj,loadfactor,loadfactor_lim,time_counter_val,Modified_NR_flag);    


else 

    disp("Check your IsProj variable")

end




end


