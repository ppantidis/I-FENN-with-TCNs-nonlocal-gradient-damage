clc
clear all
close all
format compact
 
% -------------------------------------------------------------------------
% Load Latex plotting format
func_latex_plotting_format()

% -------------------------------------------------------------------------
% Define the loadfactor value where the Local/NN solver is applied
loadfactor_lim = 0.011;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======================= INCLUDE GLOBAL VARIABLES ========================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
func_include_flags;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ============================ LOAD INPUT FILE ============================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model_name = "SNS_Coarse";
infile     = fopen(model_name + "_parameters.txt",'r');
[SolverID,TangentID,Modified_NR_flag,Staggered_flag,ncoord,ndof,quad_elem_order,lc, ... 
            increment,inc_success_counter,min_iter,max_iter,max_accept_iter, ... 
            loadfactor,dlfactor,dlfactor_incr_threshold,increment_plot_threshold,loadfactor_plot_threshold, ...
            flaglf,countflaglf,incrflag,flagplot, ...
            ndomains,nprops,materialprops,estar_shear_flag,alpha_val,beta_val,e_delta,dmax, ...
            nnodes,coords,nelem,maxnodes,connect,nelnodes,elident_vec,nfix,fixnodes] = func_read_input_file(infile);  

fclose(infile);

% -------------------------------------------------------------------------
% Find the number of elements attached to each node
num_elem_at_node = zeros(nnodes,1);
for i = 1:size(connect,1)
    for j = 1:size(connect,2)
        num_elem_at_node(connect(i,j)) = num_elem_at_node(connect(i,j)) + 1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ================== SPECIFY VALUES FOR GLOBAL VARIABLES ==================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Non-local gradient parameter (g = lc^2/2)
g = lc^2/2;

% -------------------------------------------------------------------------
% Element length
lelem = coords(1,2) - coords(1,1);

% -------------------------------------------------------------------------
% Damage approach (local vs non-local)
if SolverID == 1
    solver_string = "FEM_Local";
elseif SolverID == 2
    solver_string = "FEM_Gradient";
elseif SolverID == 3
    solver_string = "FEM_Integral";
else
    disp("Check your SolverID - FEM Main Script")
end

% -------------------------------------------------------------------------
% Analytical vs Numerical Tangent
if TangentID == 1
    tangent_string = "Analytical";
elseif TangentID == 2
    tangent_string = "Numerical";
else
    disp("Check your TangentID - FEM_Main_Script")
end

% -------------------------------------------------------------------------
% Calculate matrix with neighbouring points and their weights for the nonlocal integral method
if SolverID == 1 || SolverID == 2
    n_hood = [];
    weights = [];
elseif SolverID == 3
    [n_hood, weights] = func_nhood_gausspts(lc);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====================== SETTING UP MATRICES/VECTORS ======================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Delastic = func_Delastic;           % Calculate the elastic constitutive matrix D (3x3 for 2D problems)
fixnodes = (sortrows(fixnodes'))';  % Ensure that the fixnodes matrix lists the prescribed displacements in ascending order of nodes and dofs (as in the global system)

% Setting up: 
npoints = func_numberofintegrationpoints(quad_elem_order);   % Number of integration points   
xilist  = func_integrationpoints(quad_elem_order);           % Positions of integration points
w       = func_integrationweights(quad_elem_order);          % Weights of integration points

% Initialize vectors and matrices at zero
dofs_stored                 = zeros(ndof*nnodes,1); % Contains the values of displacements at each node
local_strain_mat_stored     = zeros(nelem,npoints);       % Contains the values of local equivalent strain at each Gauss point
nonlocal_strain_mat_stored  = zeros(nelem,npoints);       % Contains the values of nonlocal equivalent strain at each Gauss point
history_var_mat_stored      = zeros(nelem,npoints);       % Contains the values of the history variable (damage or kappa) at each Gauss point
stress_s1_mat_stored        = zeros(nelem,npoints);       % Contains the values of 1st principal stress at each Gauss point
xcoord_GP_mat               = zeros(nelem,npoints);       % Contains the x-coordinates of all Gauss points     
ycoord_GP_mat               = zeros(nelem,npoints);       % Contains the y-coordinates of all Gauss points

% -------------------------------------------------------------------------
% Create an empty tensor to store the data needed for the next TCN prediction
Data_IFENN_GP = zeros(nelem * npoints,160,5);
save Data_IFENN_GP Data_IFENN_GP 
time_counter     = zeros(160,1);
time_counter_val = 0;
IFENN_flag       = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======================= GAUSS POINTS COORDINATES ========================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the 4x4 shape function matrix (constant for all our elements)
for integ_point = 1:npoints
    xi = xilist(:,integ_point);                 % Size: 2x1
    N(:,integ_point) = func_shapefunctions(xi,quad_elem_order); % Size: 4x1
    dNdxi = func_shapefunctionderivs(xi,quad_elem_order);
end

% Find and store the x and y coordinates of all Gauss points
for lmn = 1:nelem

    lmncoord = zeros(ncoord,maxnodes);
    % Extract nodal coordinates for the current element
    for a = 1:nelnodes(lmn)
        for i = 1:ncoord
            lmncoord(i,a) = coords(i,connect(a,lmn));
        end        
    end
    
    % Transpose lmncoord (from 2x4 to 4x2)
    lmncoord_transp = lmncoord';
    % Calculate the x and y coordinates of the Gauss points
    xcoord_GP_mat(lmn,:) = N' * lmncoord_transp(:,1);
    ycoord_GP_mat(lmn,:) = N' * lmncoord_transp(:,2);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =========== JACOBIAN FOR SHAPE FUNCTION GRADIENTS (MATRIX A) ============
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the jacobian matrix for the shape function derivatives (Matrix_A)
GP_counter = 0;
Jacobian_SF_gradient = zeros(3,3,npoints*nelem); % Size: 3 x 3 x total no. GPs

for lmn = 1:nelem
    
    lmncoord = zeros(ncoord,maxnodes);
    % Extract nodal coordinates for the current element
    for a = 1:nelnodes(lmn)
        for i = 1:ncoord
            lmncoord(i,a) = coords(i,connect(a,lmn));
        end        
    end
    
    for integ_point = 1:npoints
        GP_counter  = GP_counter + 1;
        xi          = xilist(:,integ_point);                             % Size: 2x1
        dNdxi       = func_shapefunctionderivs(xi,quad_elem_order);   % Size: 4x2
        dxdxi       = lmncoord * dNdxi; % Size: 2x2 - Jacobian matrix
        
        % Note for notation
        % dxdxi = [a b; c d], where: 
        a = dxdxi(1,1); % a = partial(x)/partial(xi)
        b = dxdxi(1,2); % b = partial(y)/partial(xi)
        c = dxdxi(2,1); % c = partial(x)/partial(eta)
        d = dxdxi(2,2); % d = partial(y)/partial(eta)

        Jacobian_SF_gradient(:,:,GP_counter) = [a^2 b^2 2*a*b; c^2 d^2 2*c*d; a*c b*d a*d+c*b];

    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======================= NEWTON - RAPHSON ANALYSIS =======================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while countflaglf == 0 && incrflag < 100 && loadfactor <= 0.801
    
    % Condition to terminate the while loop: the total load is applied
    if flaglf == true; countflaglf = 1; end
    
    % Condition to terminate the while loop: dlfactor is too small
    if dlfactor < 10^-30; break; end


    % ---------------------------------------------------------------------
    % --------------------------- SOLVER ROUTINE --------------------------
    % ---------------------------------------------------------------------
      
    IsProj = 0;
    disp('---------------------------------------------------------------')
    disp('------------------------ NEW INCREMENT ------------------------')
    disp('---------------------------------------------------------------')

    disp(strcat("Loadfactor: " + loadfactor))

    % Perform the Newton Raphson analysis
    tic
    if Staggered_flag == 0
        [J, DRdu, dofs, Res_F_F_norm, Res_u_norm, last_iteration, history_var_mat, gausspoints_prop_mat, ~, Reactions_y, Res_e_norm, time_counter_val, IFENN_flag] = func_NewtonRaphson_DispControl(dofs_stored,[fixnodes(1:2,:); fixnodes(3,:)*loadfactor],increment,Delastic,history_var_mat_stored,num_elem_at_node,n_hood,weights,IsProj,model_name,loadfactor,loadfactor_lim,Modified_NR_flag);
    elseif Staggered_flag == 1
        [J, DRdu, dofs, Res_F_F_norm, Res_u_norm, last_iteration, history_var_mat, gausspoints_prop_mat, ~, Reactions_y, Res_e_norm, time_counter_val] = func_NewtonRaphson_DispControl_staggered(dofs_stored,[fixnodes(1:2,:); fixnodes(3,:)*loadfactor],increment,Delastic,history_var_mat_stored,num_elem_at_node,n_hood,weights,IsProj,model_name,loadfactor,loadfactor_lim,Modified_NR_flag);    
    else
        disp ('Check your Staggered_flag!')
    end

    if IFENN_flag == 0
        time_counter(increment,1) = toc;             
        time_counter(increment,1)                       %#ok<NOPTS> 
    elseif IFENN_flag == 1
        time_counter(increment,1) = time_counter_val; 
        time_counter(increment,1)                       %#ok<NOPTS> 
        solver_string = "IFENN";
    else 
        disp ('Check your time counter!')
    end

    % -----------------------------------------------------------------
    % ---------------- SAVING AND PLOTTING ROUTINES -------------------
    % -----------------------------------------------------------------

    if (last_iteration ~= min_iter) && (loadfactor > loadfactor_plot_threshold) % max_accept_iter instead of min_iter

        Partial_nonlocal_physical = zeros(3,1,GP_counter);

        if quad_elem_order == 2 && ndof == 3
            % -------------------------------------------------------------        
            % Compute the partial derivatives of the nonlocal strain wrt the natural coordinates (Matrix_C)
            Partial_nonlocal_natural = func_partial_nonlocal_natural(dofs);
    
            % -------------------------------------------------------------        
            % Compute the partial derivatives of the nonlocal strain wrt the physical coordinates (Matrix_B = Matrix_C ./ Matrix_A)            
            for GP = 1:GP_counter
                Partial_nonlocal_physical(:,1,GP) = inv(Jacobian_SF_gradient(:,:,GP)) * Partial_nonlocal_natural(:,1,GP);
            end
            Partial_nonlocal_physical = squeeze(Partial_nonlocal_physical);
        end

        % -----------------------------------------------------------------
        % Store reactions
        Reactions_mat(inc_success_counter,1:2) = [loadfactor, Reactions_y];

        % -----------------------------------------------------------------
        % Activate plotting mode
        IsProj = 1;
        [~, ~, ~, ~, ~, ~, ~, gausspoints_prop_mat, nodes_prop_mat, ~, ~, ~] = func_NewtonRaphson_DispControl(dofs,[fixnodes(1:2,:); fixnodes(3,:)*loadfactor],increment,Delastic,history_var_mat_stored,num_elem_at_node,n_hood,weights,IsProj,model_name,loadfactor,loadfactor_lim,Modified_NR_flag);                 

        % -------------------------------------------------------------
        % Store properties in Data_IFENN_GP
        % Save results: GAUSS POINTS 
        xcoord_GP_vec   = reshape(xcoord_GP_mat',1,[ ])';
        ycoord_GP_vec   = reshape(ycoord_GP_mat',1,[ ])';
        loadfactor_GP_vec = ones(length(xcoord_GP_vec),1) * loadfactor;
        Data_IFENN_GP(:,increment,:) = [xcoord_GP_vec ycoord_GP_vec gausspoints_prop_mat(:,8) gausspoints_prop_mat(:,9) loadfactor_GP_vec];
        save Data_IFENN_GP Data_IFENN_GP
 
%         disp(strcat('Max strain: ', num2str(max(max(history_var_mat)))))
        
        disp(strcat('Max damage: ', num2str(max(gausspoints_prop_mat(:,1)))))      

        % -----------------------------------------------------------------
        % Save results
        func_savingresults(xcoord_GP_mat,ycoord_GP_mat,g,loadfactor,gausspoints_prop_mat,increment,model_name,solver_string,tangent_string,nnodes,lelem,nodes_prop_mat,Res_F_F_norm,Res_u_norm,Data_IFENN_GP,Res_e_norm,Reactions_mat,coords,Partial_nonlocal_physical,time_counter,Modified_NR_flag)

%         % -----------------------------------------------------------------
%         % Plot at Gauss points (scatter)
%         func_plotmesh_gps(xcoord_GP_vec,ycoord_GP_vec,gausspoints_prop_mat(:,1),gausspoints_prop_mat(:,8),gausspoints_prop_mat(:,9),model_name,last_iteration,increment,loadfactor,Res_F_F_norm,Res_u_norm,solver_string,tangent_string);

%         % -----------------------------------------------------------------
%         % Plot at elements (interpolation based on nodal values)
%         func_plotmesh_elem(coords,connect,nelem,nodes_prop_mat(:,1),nodes_prop_mat(:,8),nodes_prop_mat(:,9),nelnodes,'interp',model_name,inc_success_counter,last_iteration,increment,loadfactor,Res_F_F_norm,Res_u_norm,solver_string,tangent_string);

    end
 
    % ---------------------------------------------------------------------
    % --------------------- ADAPTIVE LOADING ROUTINE ----------------------
    % ---------------------------------------------------------------------
    
    % ---------------------------------------------------------------------
    % Adapt load incrementation w.r.t. the number of iterations needed for convergence

%     if last_iteration <= min_iter                                               % Fast Convergence - dlfactor increases
%         [dofs_stored,loadfactor_stored,history_var_mat_stored,incrflag,flagplot,loadfactor,increment,inc_success_counter,flaglf,dlfactor]   = func_adapt_quick_convergence(dofs,loadfactor,history_var_mat,last_iteration,dlfactor,increment,inc_success_counter);
%     elseif (min_iter < last_iteration) && (last_iteration <= max_iter)          % Moderate Convergence - dlfactor remains the same
%         [dofs_stored,loadfactor_stored,history_var_mat_stored,incrflag,flagplot,loadfactor,increment,inc_success_counter,flaglf]            = func_adapt_moderate_convergence(dofs,loadfactor,history_var_mat,dlfactor,increment,inc_success_counter);
%     elseif (max_iter < last_iteration) && (last_iteration < max_accept_iter)    % Slow Convergence - dlfactor decreases
%         [dofs_stored,loadfactor_stored,history_var_mat_stored,incrflag,flagplot,loadfactor,increment,inc_success_counter,flaglf,dlfactor]   = func_adapt_slow_convergence(dofs,loadfactor,history_var_mat,last_iteration,dlfactor,increment,inc_success_counter);
%     else                                                                        % No Convergence - discard the last step and repeat with smaller load value
%         [flagplot,flaglf,countflaglf,incrflag,loadfactor,dlfactor]                                                                          = func_adapt_no_convergence(incrflag,loadfactor_stored,dlfactor);
%     end

    if last_iteration < max_accept_iter  % Convergence - dlfactor remains the same
        [dofs_stored,loadfactor_stored,history_var_mat_stored,incrflag,flagplot,loadfactor,increment,inc_success_counter,flaglf] = func_adapt_moderate_convergence(dofs,loadfactor,history_var_mat,dlfactor,increment,inc_success_counter);
    else                                 % No Convergence - discard the last step and repeat with smaller load value
        [flagplot,flaglf,countflaglf,incrflag,loadfactor,dlfactor] = func_adapt_no_convergence(incrflag,loadfactor_stored,dlfactor);
    end
    
end



