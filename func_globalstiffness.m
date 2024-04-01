function [J, DRdu, Res_F, history_var_mat, gausspoints_prop_mat, nodes_prop_mat, Res_e, K, time_counter_val] = func_globalstiffness(dofs,Delastic,model_name,increment,iteration,history_var_mat_previousinc,num_elem_at_node,n_hood,weights,IsProj,loadfactor,loadfactor_lim,time_counter_val,Modified_NR_flag)     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====================== ASSEMBLE THE TANGENT MATRIX ======================
% ===================== ASSEMBLE THE RESIDUAL VECTOR ======================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Include global variables
func_include_flags;

% Setting up: 
npoints = func_numberofintegrationpoints(quad_elem_order);   % Number of integration points

load Data_IFENN_GP
%%

if loadfactor < loadfactor_lim
    
    % ---------------------------------------------------------------------
    % Put a placeholder for the residual of nonlocal strain
    Res_e = NaN;
    
    % -------------------------------------------------------------------------
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOLVER ROUTINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % -------------------------------------------------------------------------
    if IsProj == 0
    
        % Create empty entries for element/nodal properties 
        gausspoints_prop_mat = [];
        nodes_prop_mat = [];

        
        % Initializing global stiffness, global internal force and gather matrices 
        % at zero:
        J                   = zeros(ndof*nnodes,ndof*nnodes);
        K                   = zeros(ndof2*nnodes,ndof2*nnodes);
        DRdu                = zeros(ndof*nnodes,ndof*nnodes);
        Res_F               = zeros(ndof*nnodes,1);
        residual_elpos      = zeros(maxnodes*ndof,maxnodes*ndof);
        residual_elneg      = zeros(maxnodes*ndof,maxnodes*ndof);
        lmncoord            = zeros(ncoord,maxnodes);
        lmndof              = zeros(ndof,maxnodes);
        history_var_mat     = zeros(nelem,npoints);
        
        % Loop over all the elements
        for lmn = 1:nelem
        
            % Extract coords of nodes, DOF for the current element
            for a = 1:nelnodes(lmn)
                for i = 1:ncoord
                    lmncoord(i,a) = coords(i,connect(a,lmn));
                end
                for i = 1:ndof
                    lmndof(i,a) = dofs(ndof*(connect(a,lmn)-1)+i);
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ANALYTICAL CALCULATION - TANGENT MATRIX (for each element): 
            
            if TangentID == 1
    
                % -------------------------------------------------------------
                if SolverID == 1 
                    [j_el, Res_F_el, history_var_mat(lmn,:), ~, ~, k_el] = func_elstif_Local(lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),alpha_val,beta_val,e_delta,dmax,n_hood,weights,1,IsProj,quad_elem_order,ndof,estar_shear_flag);
                % -------------------------------------------------------------
                elseif SolverID == 2
                    [j_el, Res_F_el, history_var_mat(lmn,:), ~, ~, k_el] = func_elstif_Nonlocgradient(lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),g,alpha_val,beta_val,e_delta,dmax,n_hood,weights,1,IsProj,quad_elem_order,ndof,ndof2,estar_shear_flag);
                end
    
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % NUMERICAL APPROXIMATION - TANGENT MATRIX (for each element):
        
            if TangentID == 2
                    
                [~, Res_F_el, history_var_mat(lmn,:), ~, ~, ~] = func_elstif_Nonlocgradient(lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),g,alpha_val,beta_val,e_delta,dmax,n_hood,weights,1,IsProj,quad_elem_order,ndof,ndof2,estar_shear_flag);
                
                % element stiffness (based on residual at plus/minus lmndof - central difference)
                dl = 1e-10; % Tolerance 
                col_id = 0; % Counter of columns in dRdu (range: 1-12)
            
                for a = 1:nelnodes(lmn)
                    for i = 1:ndof
                        
                        % Increase counter
                        col_id = col_id + 1;
            
                        % Compute residual_el at (u+dl(dof)) of size: 12x1, and store 
                        % it in the next column of residual
                        lmndofi = lmndof;    
                        lmndofi(i,a) = lmndofi(i,a) + dl;
                        
                        [~, residual_elpos(:,col_id), ~, ~, ~, ~] = func_elstif_Nonlocgradient(lmncoord,lmndofi,Delastic,history_var_mat_previousinc(lmn,:),g,alpha_val,beta_val,e_delta,dmax,n_hood,weights,0,IsProj,quad_elem_order,ndof,ndof2,estar_shear_flag); % Final size: 12x12
                        
                        % Compute residual_el at (u-dl(dof)) of size: 12x1, and store 
                        % it in the next column of residual
                        lmndofi = lmndof;    
                        lmndofi(i,a) = lmndofi(i,a) - dl;
                        
                        [~, residual_elneg(:,col_id), ~, ~, ~, ~] = func_elstif_Nonlocgradient(lmncoord,lmndofi,Delastic,history_var_mat_previousinc(lmn,:),g,alpha_val,beta_val,e_delta,dmax,n_hood,weights,0,IsProj,quad_elem_order,ndof,ndof2,estar_shear_flag); % Final size: 12x12
                        
                    end
                end
     
                % Compute "partial R over partial u" for all dofs at one step 
                dRdu_el = 1/(2*dl) * (residual_elpos - residual_elneg); % Size: 12x12
    
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ASSEMBLE GLOBAL MATRICES: 
            % a) analytical tangent K matrix
            % b) numerical tangent DRdu matrix
            % c) residual vector
            
            % a) Analytical tangent matrix
            if TangentID == 1
                for a = 1:nelnodes(lmn)
                    for i = 1:ndof
                        for b = 1:nelnodes(lmn)
                            for k = 1:ndof
                                rw = ndof*(connect(a,lmn)-1)+i;
                                cl = ndof*(connect(b,lmn)-1)+k;
                                J(rw,cl) = J(rw,cl) + j_el(ndof*(a-1)+i,ndof*(b-1)+k);
                            end
                        end
                    end
                end
            end
        
            % b) Numerical tangent matrix
            if TangentID  == 2
                for a = 1:nelnodes(lmn)
                    for i = 1:ndof
                        for b = 1:nelnodes(lmn)
                            for k = 1:ndof
                                rw = ndof*(connect(a,lmn)-1)+i;
                                cl = ndof*(connect(b,lmn)-1)+k;
                                DRdu(rw,cl) = DRdu(rw,cl) + dRdu_el(ndof*(a-1)+i,ndof*(b-1)+k);
                            end
                        end
                    end
                end
            end
        
            % c) Residual vector
            for a = 1:nelnodes(lmn)
                for i = 1:ndof
                    rw = ndof*(connect(a,lmn)-1)+i;
                    Res_F(rw,1) = Res_F(rw,1) + Res_F_el(ndof*(a-1)+i,1);
                end
            end

            % d) Assemble K
            if TangentID == 1
                for a = 1:nelnodes(lmn)
                    for i = 1:ndof2
                        for b = 1:nelnodes(lmn)
                            for k = 1:ndof2
                                rw = ndof2*(connect(a,lmn)-1)+i;
                                cl = ndof2*(connect(b,lmn)-1)+k;
                                K(rw,cl) = K(rw,cl) + k_el(ndof2*(a-1)+i,ndof2*(b-1)+k);
                            end
                        end
                    end
                end
            end
        end
    
    elseif IsProj == 1
        
    % -------------------------------------------------------------------------
    %%%%%%%%%%%%%%%%%%%%%%% PLOTTING PROPERTIES ROUTINE %%%%%%%%%%%%%%%%%%%%%%%
    % -------------------------------------------------------------------------
    
        % Create empty entries for solver variables
        J         = [];
        K         = []; 
        DRdu      = [];
        Res_F     = [];
        history_var_mat = [];
    
        % Initializing element/nodal properties and gather matrices at zero
        gausspoints_prop_mat = zeros(nelem * npoints,9);
        nodes_prop_mat       = zeros(nnodes,9);
        lmncoord             = zeros(ncoord,maxnodes);
        lmndof               = zeros(ndof,maxnodes);
    
        % Loop over all the elements
        for lmn = 1:nelem
        
            % Extract coords of nodes, DOF for the current element
            for a = 1:nelnodes(lmn)
                for i = 1:ncoord
                    lmncoord(i,a) = coords(i,connect(a,lmn));
                end
                for i = 1:ndof
                    lmndof(i,a) = dofs(ndof*(connect(a,lmn)-1)+i);
                end
            end
    
            % Compute element/nodal properties
            % -------------------------------------------------------------
            if SolverID == 1 
                [~, ~, ~, gausspoints_prop_mat(npoints*lmn-(npoints-1):npoints*lmn,:), nodes_prop_mat_elem, ~] = func_elstif_Local(lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),alpha_val,beta_val,e_delta,dmax,n_hood,weights,1,IsProj,quad_elem_order,ndof,estar_shear_flag);
            % -------------------------------------------------------------
            elseif SolverID == 2
                [~, ~, ~, gausspoints_prop_mat(npoints*lmn-(npoints-1):npoints*lmn,:), nodes_prop_mat_elem, ~] = func_elstif_Nonlocgradient(lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),g,alpha_val,beta_val,e_delta,dmax,n_hood,weights,1,IsProj,quad_elem_order,ndof,ndof2,estar_shear_flag);
            end

            % Populate the properties matrix
            for a = 1:nelnodes(lmn)
                rw = connect(a,lmn);
                nodes_prop_mat(rw,:) = nodes_prop_mat(rw,:) + nodes_prop_mat_elem(a,:);
            end
    
        end
    
        % Calculate the final matrix of nodal properties
        nodes_prop_mat = nodes_prop_mat ./ num_elem_at_node;
    
    else 
    
        disp("Check your IsProj variable - globalstiffness")
    
    end

%%
elseif loadfactor > loadfactor_lim
    % -------------------------------------------------------------------------
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRE-PROCESS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % -------------------------------------------------------------------------

    % -------------------------------------------------------------------------
    lmncoord            = zeros(ncoord,maxnodes);
    lmndof              = zeros(ndof,maxnodes);
    local_strain_mat    = zeros(nelem,npoints);
    xcoord_GP_mat       = zeros(nelem,npoints);               % Contains the x-coordinates of all Gauss points     
    ycoord_GP_mat       = zeros(nelem,npoints);               % Contains the y-coordinates of all Gauss points
    nodes_prop_mat      = zeros(nnodes,1);              % Contains the local strain at the nodes
    lelem               = coords(1,2) - coords(1,1);    % Element length
    
    % -------------------------------------------------------------------------
    % Compute the 4x4 shape function matrix (constant for all our elements)
    xilist  = func_integrationpoints(quad_elem_order);  % Positions of integration points in parent element
    for integ_point = 1:npoints
        xi = xilist(:,integ_point);                 % Size: 2x1
        N(:,integ_point) = func_shapefunctions(xi,quad_elem_order); % Size: 4x1
    end
    
    % Loop over all the elements
    for lmn = 1:nelem
    
        % Extract coords of nodes, DOF for the current element
        for a = 1:nelnodes(lmn)
            for i = 1:ncoord
                lmncoord(i,a) = coords(i,connect(a,lmn));
            end
            for i = 1:ndof
                lmndof(i,a) = dofs(ndof*(connect(a,lmn)-1)+i);
            end
        end
    
        % -----------------------------------------------------------------
        % LOCAL STRAIN VECTOR (for each element - GPs and nodes): 
        [local_strain_el, nodes_prop_vec_elem] = func_elocalvec(lmncoord,lmndof,quad_elem_order,estar_shear_flag);
        local_strain_mat(lmn,:) = local_strain_el;
        % Populate the properties matrix
        for a = 1:nelnodes(lmn)
            rw = connect(a,lmn);
            nodes_prop_mat(rw,:) = nodes_prop_mat(rw,:) + nodes_prop_vec_elem(a,:);
        end

        % -----------------------------------------------------------------
        % Transpose lmncoord (from 2x4 to 4x2) and calculate the x and y 
        % coordinates of the Gauss points
        lmncoord_transp = lmncoord';
        xcoord_GP_mat(lmn,:) = N' * lmncoord_transp(:,1);
        ycoord_GP_mat(lmn,:) = N' * lmncoord_transp(:,2);
    
    end
    
    % ---------------------------------------------------------------------
    % GAUSS POINTS PROPERTIES MATRIX 
    xcoord_GP_vec               = reshape(xcoord_GP_mat',1,[ ])';
    ycoord_GP_vec               = reshape(ycoord_GP_mat',1,[ ])';
    local_strain_vec            = reshape(local_strain_mat',1,[ ])';
    loadfactor_GP_vec           = ones(length(xcoord_GP_vec),1) * loadfactor;

    data_StrB100_GP        = zeros(nelem * npoints,1,5);
    data_StrB100_GP(:,1,:) = [xcoord_GP_vec ycoord_GP_vec local_strain_vec zeros(length(xcoord_GP_vec),1) loadfactor_GP_vec];    
    
    Data_IFENN_GP(:,increment,:) = data_StrB100_GP;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    time_counter_val = toc + time_counter_val;
    matfilename            = sprintf('PredData_IFENN_GP_time_features_lf_%d_%d.mat', increment, iteration);
    save(matfilename, 'Data_IFENN_GP'); 

    % ---------------------------------------------------------------------
    % Execute Python to predict the non-local strain and all the derivatives
    pyrunfile("./TCN_IFENN.py", increment = int2str(increment), iteration = int2str(iteration));

    % Load the predictions in MATLAB
    load IFENN_Predictions_SNS_Coarse_D10_P00_dec10000_Ad2000_LBFGS5000.mat
    
    % ---------------------------------------------------------------------
    % Slice the increment of interest and assign to vectors
    % For Full Newton Raphson (predictions and derivatives)l
    if Modified_NR_flag == 0
        enonlocal_strain_vec      = predictions(:,increment); 
        denonlocal_delocal_vec    = enonlocal_elocal(:,increment);
%         denonlocal_gradient_vec   = squeeze(enonlocal_grad_xy(:,increment,:));
        denonlocal_dx_vec         = zeros(nelem * npoints, 1); 
        denonlocal_dy_vec         = zeros(nelem * npoints, 1);   
        
%         % -----------------------------------------------------------------
%         matfilename1 = sprintf('enonlocal_strain_vec_%d_%d.mat', increment, iteration);
%         matfilename2 = sprintf('denonlocal_delocal_vec_%d_%d.mat', increment, iteration);
%         matfilename3 = sprintf('denonlocal_dx_vec_%d_%d.mat', increment, iteration);
%         matfilename4 = sprintf('denonlocal_dy_vec_%d_%d.mat', increment, iteration);
%         save(matfilename1, 'enonlocal_strain_vec');    
%         save(matfilename2, 'denonlocal_delocal_vec');
%         save(matfilename3, 'denonlocal_dx_vec'); 
%         save(matfilename4, 'denonlocal_dy_vec'); 
        
    % For Modified Newton Raphson (predictions only)
    elseif Modified_NR_flag == 1
        enonlocal_strain_vec      = predictions(:,increment); 
        denonlocal_delocal_vec    = zeros(nelem * npoints, 1);
        denonlocal_dx_vec         = zeros(nelem * npoints, 1);
        denonlocal_dy_vec         = zeros(nelem * npoints, 1);   
        
%         % -----------------------------------------------------------------
%         matfilename1 = sprintf('enonlocal_strain_vec_%d_%d.mat', increment, iteration);
%         save(matfilename1, 'enonlocal_strain_vec');    
        
    end
    tic 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % -------------------------------------------------------------------------
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOLVER ROUTINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % -------------------------------------------------------------------------
    if IsProj == 0
    
        % Create empty entries for element/nodal properties 
        gausspoints_prop_mat = [];
        nodes_prop_mat = [];
       
        % Initializing global stiffness, global internal force and gather matrices 
        % at zero:
        J                   = zeros(ndof*nnodes,ndof*nnodes);
        K                   = zeros(ndof2*nnodes,ndof2*nnodes);
        DRdu                = zeros(ndof*nnodes,ndof*nnodes);
        Res_F               = zeros(ndof*nnodes,1);
        Res_e               = zeros(nnodes,1);        
        history_var_mat     = zeros(nelem,npoints);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ---------------------------------------------------------------------
        % Loop over all the elements
        for lmn = 1:nelem
        
            % Extract coords of nodes, DOF for the current element
            for a = 1:nelnodes(lmn)
                for i = 1:ncoord
                    lmncoord(i,a) = coords(i,connect(a,lmn));
                end
                for i = 1:ndof
                    lmndof(i,a) = dofs(ndof*(connect(a,lmn)-1)+i);
                end
            end
            
            % Extract the nonlocal strain and the derivative for the specific element
            enonlocal_strain_el     = enonlocal_strain_vec(npoints*lmn-(npoints-1):npoints*lmn)';
            denonlocal_delocal_el   = denonlocal_delocal_vec(npoints*lmn-(npoints-1):npoints*lmn)';
            denonlocal_dx_el        = denonlocal_dx_vec(npoints*lmn-(npoints-1):npoints*lmn)';
            denonlocal_dy_el        = denonlocal_dy_vec(npoints*lmn-(npoints-1):npoints*lmn)';
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ANALYTICAL CALCULATION - TANGENT MATRIX (for each element): 
            
            if TangentID == 1
%                 [j_el, Res_F_el, Res_e_el, history_var_mat(lmn,:), ~, ~, k_el] = func_elstif_Local_NN_history_kappa(lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),g,alpha_val,beta_val,e_delta,dmax,1,IsProj,enonlocal_strain_el,denonlocal_delocal_el,denonlocal_dx_el,denonlocal_dy_el,quad_elem_order,ndof,Modified_NR_flag,estar_shear_flag);
                [j_el, Res_F_el, Res_e_el, history_var_mat(lmn,:), ~, ~, k_el] = func_elstif_Local_NN(lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),g,alpha_val,beta_val,e_delta,dmax,1,IsProj,enonlocal_strain_el,denonlocal_delocal_el,denonlocal_dx_el,denonlocal_dy_el,quad_elem_order,ndof,Modified_NR_flag,estar_shear_flag);
            else
                disp("Check your TangentID - globalstiffness")
            end
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ASSEMBLE GLOBAL MATRICES: 
            % a) analytical tangent K matrix
            % b) residual vector
    
            % a) Analytical tangent matrix
            for a = 1:nelnodes(lmn)
                for i = 1:ndof
                    for b = 1:nelnodes(lmn)
                        for k = 1:ndof
                            rw = ndof*(connect(a,lmn)-1)+i;
                            cl = ndof*(connect(b,lmn)-1)+k;
                            J(rw,cl) = J(rw,cl) + j_el(ndof*(a-1)+i,ndof*(b-1)+k);
                        end
                    end
                end
            end
        
            % b) Residual vector
            for a = 1:nelnodes(lmn)
                for i = 1:ndof
                    rw = ndof*(connect(a,lmn)-1)+i;
                    Res_F(rw,1) = Res_F(rw,1) + Res_F_el(ndof*(a-1)+i,1);
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%% TRIPLE CHECK THE COMPUTATION BELOW %%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % c) Nonlocal strain PDE residual vector
            for a = 1:nelnodes(lmn)
                rw = (connect(a,lmn)-1)+ndof1;
                Res_e(rw,1) = Res_e(rw,1) + Res_e_el(ndof1*(a-1)+1,1);
            end

            % d) Analytical stiffness matrix
            for a = 1:nelnodes(lmn)
                for i = 1:ndof
                    for b = 1:nelnodes(lmn)
                        for k = 1:ndof
                            rw = ndof*(connect(a,lmn)-1)+i;
                            cl = ndof*(connect(b,lmn)-1)+k;
                            K(rw,cl) = K(rw,cl) + k_el(ndof*(a-1)+i,ndof*(b-1)+k);
                        end
                    end
                end
            end

        end
    
    
    elseif IsProj == 1
    % -------------------------------------------------------------------------
    %%%%%%%%%%%%%%%%%%%%%%% PLOTTING PROPERTIES ROUTINE %%%%%%%%%%%%%%%%%%%%%%%
    % -------------------------------------------------------------------------
    
        % Create empty entries for solver variables
        J         = [];
        K         = [];     
        DRdu      = [];
        Res_F     = [];
        Res_e     = [];        
        history_var_mat = [];
    
        % Initializing element/nodal properties and gather matrices at zero
        gausspoints_prop_mat = zeros(nelem * func_numberofintegrationpoints(quad_elem_order),9);
        nodes_prop_mat       = zeros(nnodes,9);
        lmncoord             = zeros(ncoord,maxnodes);
        lmndof               = zeros(ndof,maxnodes);
    
        % Loop over all the elements
        for lmn = 1:nelem
        
            % Extract coords of nodes, DOF for the current element
            for a = 1:nelnodes(lmn)
                for i = 1:ncoord
                    lmncoord(i,a) = coords(i,connect(a,lmn));
                end
                for i = 1:ndof
                    lmndof(i,a) = dofs(ndof*(connect(a,lmn)-1)+i);
                end
            end
    
            % Extract local strain, nonlocal strain and derivative for the
            % specific element
            enonlocal_strain_el    = enonlocal_strain_vec(npoints*lmn-(npoints-1):npoints*lmn)';
            denonlocal_delocal_el  = denonlocal_delocal_vec(npoints*lmn-(npoints-1):npoints*lmn)';
            denonlocal_dx_el       = denonlocal_dx_vec(npoints*lmn-(npoints-1):npoints*lmn)';
            denonlocal_dy_el       = denonlocal_dy_vec(npoints*lmn-(npoints-1):npoints*lmn)';

            % Compute element/nodal properties
%             [~, ~, ~, ~, gausspoints_prop_mat(npoints*lmn-(npoints-1):npoints*lmn,:), nodes_prop_mat_elem, ~] = func_elstif_Local_NN_history_kappa(lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),g,alpha_val,beta_val,e_delta,dmax,1,IsProj,enonlocal_strain_el,denonlocal_delocal_el,denonlocal_dx_el,denonlocal_dy_el,quad_elem_order,ndof,Modified_NR_flag,estar_shear_flag);  
            [~, ~, ~, ~, gausspoints_prop_mat(npoints*lmn-(npoints-1):npoints*lmn,:), nodes_prop_mat_elem, ~] = func_elstif_Local_NN(lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),g,alpha_val,beta_val,e_delta,dmax,1,IsProj,enonlocal_strain_el,denonlocal_delocal_el,denonlocal_dx_el,denonlocal_dy_el,quad_elem_order,ndof,Modified_NR_flag,estar_shear_flag);              
            
            % Populate the properties matrix
            for a = 1:nelnodes(lmn)
                rw = connect(a,lmn);
                nodes_prop_mat(rw,:) = nodes_prop_mat(rw,:) + nodes_prop_mat_elem(a,:);
            end
    
        end
    
        % Calculate the final matrix of nodal properties
        nodes_prop_mat = nodes_prop_mat ./ num_elem_at_node;
    
    else 
    
        disp("Check your IsProj variable - globalstiffness")
    
    end


else

    disp("Check your loadfactor variable - globalstiffness")


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end