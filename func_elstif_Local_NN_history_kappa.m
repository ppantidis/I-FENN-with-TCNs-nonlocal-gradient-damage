function [j_el, f_internal_el, res_e_el, kappa_el, gausspoints_prop_mat_elem, nodes_prop_mat_elem, k_el] = func_elstif_Local_NN_history_kappa(lmncoord,dofs,Delastic,kappa_el_previousinc,g,alpha_val,beta_val,e_delta,dmax,IsM,IsProj,enonlocal_strain_el,denonlocal_delocal_el,denonlocal_dx_el,denonlocal_dy_el,quad_elem_order,ndof,Modified_NR_flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ============ ELEMENT STIFFNESS MATRIX & INTERNAL FORCE VECTOR ===========
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

denonlocal_gradient_el = [denonlocal_dx_el; denonlocal_dy_el];

if quad_elem_order == 1
    number_of_nodes = 4;
elseif quad_elem_order == 2
    number_of_nodes = 8;
else
    disp('Check your quad_elem_order variable')
end

% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOLVER ROUTINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
if IsProj == 0

    % Create empty entries for element/nodal properties 
    gausspoints_prop_mat_elem = [];
    nodes_prop_mat_elem = [];   

    % Isolate x-y dofs 
    dofs_u = dofs(1:2,:); % 8x1
    
    % Unroll dofs into a vector
    dofs_vec = dofs_u(:); % 8x1

    % Setting up: 
    npoints = func_numberofintegrationpoints(quad_elem_order);   % Number of integration points
    xilist  = func_integrationpoints(quad_elem_order);           % Positions of integration points
    w       = func_integrationweights(quad_elem_order);          % Weights of integration points
    
    % Initialize matrices and vectors at zero
    j_el            = zeros(number_of_nodes*ndof,number_of_nodes*ndof);
    k_el            = zeros(number_of_nodes*ndof,number_of_nodes*ndof);
    f_internal_el   = zeros(number_of_nodes*ndof,1);
    res_e_el        = zeros(number_of_nodes,1);
    damage_el       = zeros(1,npoints);
    kappa_el        = zeros(1,npoints);

    % Loop over the integration points
    for integ_point = 1:npoints
    
        % ================== SHAPE FUNCTIONS AND DERIVATIVES ==================
        % Compute shape functions && derivatives wrt local coords
        xi    = xilist(:,integ_point);          % Size: 2x1
        N     = func_shapefunctions(xi,quad_elem_order);        % Size: 4x1 
        dNdxi = func_shapefunctionderivs(xi,quad_elem_order);   % Size: 4x2
        
        % Compute the jacobian matrix && its determinant
        dxdxi = lmncoord * dNdxi; % Size: 2x2 - Jacobian matrix
        dt    = det(dxdxi);
    
        % Convert shape function derivatives to derivatives wrt global coords
        dNdx = dNdxi/dxdxi; % Size: 4x2 - Multiplies by inverse of dxdxi
    
        % Transpose and expand dNdx (entries of B matrix) to B matrix
        % dNdx is 4x2 ===> B is 3x8
        [B, Be] = func_Bmatrix(dNdx,quad_elem_order);
    
    
        % ==================== LOCAL STRAIN CALCULATIONS ======================
        % Compute the strain by multiplying the shape function derivatives with
        % the displacements
        strain_vec_gxy = B * dofs_vec; % Size: 3x1 (3x8 * 8x1)
    
        % Extract the strain values for each integration point
        exx = strain_vec_gxy(1,1);
        eyy = strain_vec_gxy(2,1);
        gxy = strain_vec_gxy(3,1);
        
        % Calculate the equivalent strain e_star
        [e_star,s] = func_estar(exx, eyy, gxy);
        
        % Add NaN condition for [s]
        if isnan(s); s(:) = 0; end
    
        
        % ================ NONLOCAL STRAIN CALCULATIONS ===================
        enonlocal_strain    = enonlocal_strain_el(1,integ_point);
        denonlocal_delocal  = denonlocal_delocal_el(1,integ_point);
        denonlocal_gradient = denonlocal_gradient_el(:,integ_point);

        % Store the kappa variable for each Gauss Point
        % (IsM: 1 for analytical, 0 for numerical)
        if IsM == 0
            kappa_el(1,integ_point) = enonlocal_strain; 
        else
            kappa_el(1,integ_point) = max(kappa_el_previousinc(1,integ_point),enonlocal_strain);
        end


        % =============== MAZAR DAMAGE MODEL IMPLEMENTATION ===============
        % Calculate damage variable omega
        [omega, domega_denonlocal] = func_mazarmodel_Local_NN(enonlocal_strain,alpha_val,beta_val,e_delta,dmax);
        omega_previous_inc         = func_mazarmodel_Local_NN(kappa_el_previousinc(1,integ_point),alpha_val,beta_val,e_delta,dmax);
    
        % Store the damage variable for each Gauss Point
        if IsM == 0
            damage_el(1,integ_point) = omega; 
        else
            damage_el(1,integ_point) = max(omega_previous_inc,omega);
        end
        
        % ======================= STRESS CALCULATION ==========================
        stress_vec = (1 - damage_el(1,integ_point)) * Delastic * strain_vec_gxy; % Size: 3x1 (1x1 * 3x3 * 3x1)
    
        
        % ==================== ELEMENT STIFFNESS MATRIX =======================
        % Compute the contribution of this Gauss point to the element stiffness matrix
        KA = (1 - damage_el(1,integ_point)) * B' * Delastic * B;
        
        if Modified_NR_flag == 0
            KB = zeros(number_of_nodes*ndof,number_of_nodes*ndof);
            for i = 1:number_of_nodes
                for j = 1:2
                   row = 2*(i-1)+j;           
                   for ii = 1:number_of_nodes
                       for jj = 1:2
                           col = 2 * (ii - 1) + jj;
                           dD_duj = -1 * domega_denonlocal * denonlocal_delocal * (s(1) * B(1,col) + s(2) * B(2,col) + s(3) * B(3,col));     % dD_duj = - domega_dkappa * denonlocal_delocal * destar/deij * deij/duj
                           for ki = 1:number_of_nodes
                               for kj = 1:2
                                   kind = 2 * (ki - 1) + kj;
                                   KB(row,col) = KB(row,col) + dD_duj * (1 / (1 - damage_el(1,integ_point))) * KA(row,kind) * dofs_vec(kind);
                               end
                           end
                       end
                   end           
                end        
            end
        end

        % ======== Compute the contribution in the tangent matrix =========
        if Modified_NR_flag == 0
            Kint = KA + KB;
        elseif Modified_NR_flag == 1
            Kint = KA;
        else
            disp("Check your Modified Newton Raphson flag!")
        end

        % ====== Assemble the element tangent and stiffness matrices ======
        j_el  = j_el + (Kint * w(integ_point) * dt) ; % Size: 8x8 (1x1 * 8x3 * 3x3 * 3x8 * 1x1 * 1x1)
        k_el  = k_el + (KA * w(integ_point) * dt);

        % ================= ELEMENT INTERNAL FORCE VECTOR ==================
        % Compute the contribution of each Gauss Point to the element internal force
        f_internal_el = f_internal_el + B' * stress_vec * w(integ_point) * dt; % Size: 8x1

        % ================= ELEMENT PDE RESIDUAL VECTOR ===================
        % Compute the contribution of each Gauss Point to the element nonlocal strain PDE residual 
        res_e_el = res_e_el + (N * enonlocal_strain + Be' * g * denonlocal_gradient - N * e_star) * w(integ_point) * dt; % Size: 4x1
    
    end



elseif IsProj == 1
% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%% PLOTTING PROPERTIES ROUTINE %%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
    
j_el          = [];
k_el          = [];
f_internal_el = [];
res_e_el      = [];

    % Unroll dofs into a vector
    dofs_vec = dofs(:); % 8x1

    % Setting up: 
    npoints = func_numberofintegrationpoints(quad_elem_order);   % Number of integration points
    xilist  = func_integrationpoints(quad_elem_order);           % Positions of integration points
    w       = func_integrationweights(quad_elem_order);          % Weights of integration points

    % Initialize matrices and vectors at zero
    local_strain_el           = zeros(1,npoints);
    damage_el                 = zeros(1,npoints);
    kappa_el                  = zeros(1,npoints);
    M                         = zeros(number_of_nodes,number_of_nodes); 
    gausspoints_prop_mat_elem = zeros(number_of_nodes,9);
    prop_gpts                 = zeros(number_of_nodes,9);
    
    % Loop over the integration points
    for integ_point = 1:npoints
    
        % ================== SHAPE FUNCTIONS AND DERIVATIVES ==================
        % Compute shape functions && derivatives wrt local coords
        xi = xilist(:,integ_point);             % Size: 2x1
        N = func_shapefunctions(xi,quad_elem_order);            % Size: 4x1 
        dNdxi = func_shapefunctionderivs(xi,quad_elem_order);   % Size: 4x2
        
        % Compute the jacobian matrix && its determinant
        dxdxi = lmncoord * dNdxi; % Size: 2x2 - Jacobian matrix
        dt = det(dxdxi);
    
        % Convert shape function derivatives to derivatives wrt global coords
        dNdx = dNdxi/dxdxi; % Size: 4x2 - Multiplies by inverse of dxdxi
    
        % Transpose and expand dNdx (entries of B matrix) to B matrix
        % dNdx is 4x2 ===> B is 3x8
        [B, ~] = func_Bmatrix(dNdx,quad_elem_order);
    
    
        % ==================== LOCAL STRAIN CALCULATIONS ======================
        % Compute the strain by multiplying the shape function derivatives with
        % the displacements
        strain_vec_gxy = B * dofs_vec; % Size: 3x1 (3x8 * 8x1)
    
        % Extract the strain values for each integration point
        exx = strain_vec_gxy(1,1);
        eyy = strain_vec_gxy(2,1);
        gxy = strain_vec_gxy(3,1);
        exy = gxy/2;

        % Calculate the equivalent strain e_star
        [e_star,~] = func_estar(exx, eyy, gxy);
        local_strain_el(1,integ_point) = e_star;


        % ================ NONLOCAL STRAIN CALCULATIONS ===================
        enonlocal_strain    = enonlocal_strain_el(1,integ_point);

        % Store the kappa variable for each Gauss Point
        if IsM == 0
            kappa_el(1,integ_point) = enonlocal_strain; 
        else
            kappa_el(1,integ_point) = max(kappa_el_previousinc(1,integ_point),enonlocal_strain);
        end
        

        % =============== MAZAR DAMAGE MODEL IMPLEMENTATION ===============
        % Calculate damage variable omega
        [omega, ~]          = func_mazarmodel_Local_NN(enonlocal_strain,alpha_val,beta_val,e_delta,dmax);
        omega_previous_inc  = func_mazarmodel_Local_NN(kappa_el_previousinc(1,integ_point),alpha_val,beta_val,e_delta,dmax);
    
        % Store the damage variable for each Gauss Point
        if IsM == 0
            damage_el(1,integ_point) = omega; 
        else
            damage_el(1,integ_point) = max(omega_previous_inc,omega);
        end
        
        % ======================= STRESS CALCULATION ==========================
        stress_vec = (1 - damage_el(1,integ_point)) * Delastic * strain_vec_gxy; % Size: 3x1 (1x1 * 3x3 * 3x1)
    
        % -----------------------------------------------------------------
        % ========== ASSEMBLE PROPERTIES AT ELEMENT GAUSS POINT ===========
        % -----------------------------------------------------------------
        gausspoints_prop_mat_elem(integ_point,:) = [damage_el(1,integ_point) stress_vec' strain_vec_gxy' e_star enonlocal_strain];


        % -----------------------------------------------------------------
        % ===== ASSEMBLE GAUSS POINT PROPERTIES FOR NODAL PROJECTION ======
        % -----------------------------------------------------------------
        for i = 1:number_of_nodes % number of nodes per element: 4
            prop_gpts(i,1)   = prop_gpts(i,1)   + N(i) * w(integ_point) * dt * damage_el(1,integ_point); % damage
            prop_gpts(i,2:4) = prop_gpts(i,2:4) + N(i) * w(integ_point) * dt * stress_vec';              % stresses [sigma_xx,sigma_yy,tau_xy]
            prop_gpts(i,5:7) = prop_gpts(i,5:7) + N(i) * w(integ_point) * dt * [exx,eyy,exy];            % strains [e_xx, e_yy, gxy]
            prop_gpts(i,8)   = prop_gpts(i,8)   + N(i) * w(integ_point) * dt * e_star;                   % local equivalent strain [e_star]
            prop_gpts(i,9)   = prop_gpts(i,9)   + N(i) * w(integ_point) * dt * enonlocal_strain;         % nonlocal equivalent strain             
        end
    
        % Calculate M matrix for projections
        M = M + (N * N' * w(integ_point) * dt); 

    end

    % ====================== COMPUTE NODAL PROPERTIES =====================
    nodes_prop_mat_elem = M \ prop_gpts; 

else

    disp("Check your IsProj variable")

end



end

