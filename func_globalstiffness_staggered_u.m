function [J, DRdu, Res_F, history_var_mat, gausspoints_prop_mat, nodes_prop_mat, Res_e, K, time_counter_val] = func_globalstiffness_staggered_u(dofs,Delastic,model_name,increment,iteration,history_var_mat_previousinc,num_elem_at_node,n_hood,weights,IsProj,loadfactor,loadfactor_lim,time_counter_val,Modified_NR_flag)     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====================== ASSEMBLE THE TANGENT MATRIX ======================
% ===================== ASSEMBLE THE RESIDUAL VECTOR ======================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Include global variables
func_include_flags;

% Setting up: 
npoints = func_numberofintegrationpoints(quad_elem_order);   % Number of integration points
    
% ---------------------------------------------------------------------
% Put a placeholder for the residual of nonlocal strain
Res_e = NaN;


% Create empty entries for element/nodal properties 
gausspoints_prop_mat = [];
nodes_prop_mat = [];

% Initializing global stiffness, global internal force and gather matrices 
% at zero:
J                   = zeros(ndof2*nnodes,ndof2*nnodes);
K                   = zeros(ndof2*nnodes,ndof2*nnodes);
DRdu                = zeros(ndof2*nnodes,ndof2*nnodes);
Res_F               = zeros(ndof2*nnodes,1);
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

        if SolverID == 2
            [j_el, Res_F_el, history_var_mat(lmn,:), ~, ~, k_el] = func_elstif_staggered_u(lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),g,alpha_val,beta_val,e_delta,dmax,n_hood,weights,1,IsProj,quad_elem_order,ndof,ndof2,estar_shear_flag);
        else
            disp('Check the flags! - func_globalstiffness_Staggered')
            break
        end

    end
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ASSEMBLE GLOBAL MATRICES: 
    % a) analytical tangent K matrix
    % b) numerical tangent DRdu matrix
    % c) residual vector
    
    % a) Analytical tangent matrix
    if TangentID == 1
        for a = 1:nelnodes(lmn)
            for i = 1:ndof2
                for b = 1:nelnodes(lmn)
                    for k = 1:ndof2
                        rw = ndof2*(connect(a,lmn)-1)+i;
                        cl = ndof2*(connect(b,lmn)-1)+k;
                        J(rw,cl) = J(rw,cl) + j_el(ndof2*(a-1)+i,ndof2*(b-1)+k);
                    end
                end
            end
        end
    end

    % c) Residual vector
    for a = 1:nelnodes(lmn)
        for i = 1:ndof2
            rw = ndof2*(connect(a,lmn)-1)+i;
            Res_F(rw,1) = Res_F(rw,1) + Res_F_el(ndof2*(a-1)+i,1);
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


end