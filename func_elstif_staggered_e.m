function [j_el, residual_el, kappa_el, gausspoints_prop_mat_elem, nodes_prop_mat_elem] = func_elstif_staggered_e(lmncoord,dofs,Delastic,kappa_el_previousinc,g,alpha_val,beta_val,e_delta,dmax,n_hood,weights,IsM,IsProj,quad_elem_order,ndof,ndof2,estar_shear_flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ============ ELEMENT STIFFNESS MATRIX & INTERNAL FORCE VECTOR ===========
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        
    % Isolate the nodal displacements ux and uy and unroll into a vector
    dofs_u      = dofs(1:2,:);   % Size: 2x4
    dofs_u_vec  = dofs_u(:); % Size: 8x1
    
    % Isolate the nodal nonlocal equivalent strains and unroll into a vector
    dofs_e      = dofs(3,:);     % Size: 1x4
    dofs_e_vec  = dofs_e(:); % Size: 4x1

    % Setting up: 
    npoints = func_numberofintegrationpoints(quad_elem_order);   % Number of integration points
    xilist  = func_integrationpoints(quad_elem_order);           % Positions of integration points
    w       = func_integrationweights(quad_elem_order);          % Weights of integration points
    
    % Initialize matrices and vectors at zero
    kappa_el            = zeros(1,npoints);
    f_external_ee_el    = zeros(number_of_nodes,1);
    j_ee_el = 0;
    
    
    % Loop over the integration points
    for integ_point = 1:npoints
    
        % ================== SHAPE FUNCTIONS AND DERIVATIVES ==================
        % Compute shape functions && derivatives wrt local coords. We use the same shape functions for displacement and strain calculations
        xi = xilist(:,integ_point);                             % Size: 2x1
        N = func_shapefunctions(xi,quad_elem_order);            % Size: 4x1 
        dNdxi = func_shapefunctionderivs(xi,quad_elem_order);   % Size: 4x2
    
        % Compute the jacobian matrix && its determinant
        dxdxi = lmncoord * dNdxi; % Size: 2x2 - Jacobian matrix
        dt = det(dxdxi);
    
        % Convert shape function derivatives to derivatives wrt global coords
        dNdx = dNdxi/dxdxi; % Size: 4x2 
    
        % Construct B matrices from dNdx. Bu (3x8) is used for the displacements 
        % and Be (2x4) is used for the nonlocal equivalent strains
        [Bu, Be] = func_Bmatrix(dNdx,quad_elem_order);
    
     
        % ==================== LOCAL STRAIN CALCULATIONS ======================
        % Compute the local strain at the Gauss points
        strain_vec_gxy = Bu * dofs_u_vec; % Size: 3x1
        
        % Extract the strain values at each Gauss point
        exx = strain_vec_gxy(1,1);
        eyy = strain_vec_gxy(2,1);
        gxy = strain_vec_gxy(3,1);
        
        % Calculate the equivalent strain e_star and the derivatives {s}
        if estar_shear_flag == 0
            [e_star, ~] = func_estar_tension(exx, eyy, gxy);
        elseif estar_shear_flag == 1
            [e_star, ~] = func_estar_shear(exx, eyy, gxy);
        else
            disp('Check your estar_shear_flag! - func_elstif_Nonlocgradient')
        end
        
        
        % ======================== STIFFNESS MATRICES =========================
        j_ee_el = j_ee_el + (N * N' + g * (Be' * Be)) * w(integ_point) * dt;                              % Size: 4x4 (4x1 * 1x4 + 1x1 * 4x2 * 2x4 * 1x1 * 1x1)        
        j_el    = j_ee_el;

        % ============ ELEMENT NONLOCAL EXTERNAL FORCE VECTOR =================
        % Compute the contribution of each Gauss Point to the element nonlocal external force vector
        f_external_ee_el = f_external_ee_el + N * e_star * w(integ_point) * dt; % Size: 4x1 (4x1 * 1x1 * 1x1 * 1x1)
    
    end    
    
    % =============== ELEMENT INTERNAL FORCE VECTOR continued =================
    % Assemble the element internal force 
    f_internal_ee_el = j_ee_el * dofs_e_vec;              % Size: 4x1 (4x4 * 4x1)
    f_internal_el = f_internal_ee_el; % Size: 12x1
    
    
    % =============== ELEMENT EXTERNAL FORCE VECTOR continued =================
    % Assemble the element external force vector. NOTE: The zeros correspond to 
    % f_external_uu, which is the reactions values and they are calculated in 
    % the Newton-Raphson script. With this approach, we contribute with zero 
    % to those values while maintaining the same partition for k_el, 
    % f_internal, f_external. This will allow us to assemble correctly the
    % global external force vector in the globalstiffness script . 
    f_external_el = f_external_ee_el; 
    
    
    % ==================== ELEMENT RESIDUAL FORCE VECTOR ======================
    % Assemble the element residual vector.
    residual_el = [f_internal_el - f_external_el];
    
elseif IsProj == 1
% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%% PLOTTING PROPERTIES ROUTINE %%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------

    j_el = [];
    residual_el = [];
    kappa_el = [];

    % Isolate the nodal displacements ux and uy and unroll into a vector
    dofs_u      = dofs(1:2,:);   % Size: 2x4
    dofs_u_vec  = dofs_u(:); % Size: 8x1
    
    % Isolate the nodal nonlocal equivalent strains and unroll into a vector
    dofs_e      = dofs(3,:);     % Size: 1x4
    dofs_e_vec  = dofs_e(:); % Size: 4x1

    % Setting up: 
    npoints = func_numberofintegrationpoints(quad_elem_order);   % Number of integration points
    xilist  = func_integrationpoints(quad_elem_order);           % Positions of integration points
    w       = func_integrationweights(quad_elem_order);          % Weights of integration points
    
    % Initialize matrices and vectors at zero
    local_strain_el           = zeros(1,npoints);
    nonlocal_strain_el        = zeros(1,npoints);
    damage_el                 = zeros(1,npoints);
    M                         = zeros(number_of_nodes,number_of_nodes); 
    gausspoints_prop_mat_elem = zeros(number_of_nodes,9);
    prop_gpts                 = zeros(number_of_nodes,9);
    
    % Loop over the integration points
    for integ_point = 1:npoints
    
        % ================== SHAPE FUNCTIONS AND DERIVATIVES ==================
        % Compute shape functions && derivatives wrt local coords. We use the same shape functions for displacement and strain calculations
        xi = xilist(:,integ_point);             % Size: 2x1
        N = func_shapefunctions(xi,quad_elem_order);            % Size: 4x1 
        dNdxi = func_shapefunctionderivs(xi,quad_elem_order);   % Size: 4x2
    
        % Compute the jacobian matrix && its determinant
        dxdxi = lmncoord * dNdxi; % Size: 2x2 - Jacobian matrix
        dt = det(dxdxi);

        % Convert shape function derivatives to derivatives wrt global coords
        dNdx = dNdxi/dxdxi; % Size: 4x2 
    
        % Construct B matrices from dNdx. Bu (3x8) is used for the displacements 
        % and Be (2x4) is used for the nonlocal equivalent strains
        [Bu, ~] = func_Bmatrix(dNdx,quad_elem_order);
    
     
        % ==================== LOCAL STRAIN CALCULATIONS ======================
        % Compute the local strain at the Gauss points
        strain_vec_gxy = Bu * dofs_u_vec; % Size: 3x1
        
        % Extract the strain values at each Gauss point
        exx = strain_vec_gxy(1,1);
        eyy = strain_vec_gxy(2,1);
        gxy = strain_vec_gxy(3,1);
        exy = gxy/2;
        
        % Calculate the equivalent strain e_star and the derivatives {s}
        if estar_shear_flag == 0
            [e_star, ~] = func_estar_tension(exx, eyy, gxy);
        elseif estar_shear_flag == 1
            [e_star, ~] = func_estar_shear(exx, eyy, gxy);
        else
            disp('Check your estar_shear_flag! - func_elstif_Nonlocgradient')
        end
        local_strain_el(1,integ_point) = e_star;
    
        
        % ================== NON-LOCAL STRAIN CALCULATIONS ====================
        % Compute and store the nonlocal equivalent strain for each Gauss Point
        nonlocal_strain_el(1,integ_point) = N' * dofs_e_vec; % Size: 1x1
    
    
        % =================== MAZAR MODEL IMPLEMENTATION ======================
        % Calculate the following variables for each Gauss point:
        % a) nonlocal equivalent strain history parameter (kappa)
        % b) damage variable (omega) 
        % c) gradient of omega w.r.t kappa (domega_dkappa)
        [~, omega, ~] = func_mazarmodel_Nonlocgradient(nonlocal_strain_el(1,integ_point),kappa_el_previousinc(1,integ_point),alpha_val,beta_val,e_delta,dmax,IsM);
        
        % Store the output of the Mazar model to the appropriate element vectors
        damage_el(1,integ_point) = omega;
    
        
        % ======================= STRESS CALCULATION ==========================
        stress_vec = (1 - damage_el(1,integ_point)) * Delastic * strain_vec_gxy; % Size: 3x1 (1x1 * 3x3 * 3x1) 
        

        % -----------------------------------------------------------------
        % ========== ASSEMBLE PROPERTIES AT ELEMENT GAUSS POINT ===========
        % -----------------------------------------------------------------
        gausspoints_prop_mat_elem(integ_point,:) = [damage_el(1,integ_point) stress_vec' strain_vec_gxy' e_star nonlocal_strain_el(1,integ_point)];


        % -----------------------------------------------------------------
        % ===== ASSEMBLE GAUSS POINT PROPERTIES FOR NODAL PROJECTION ======
        % -----------------------------------------------------------------
        for i = 1:4 % number of nodes per element: 4
            prop_gpts(i,1)   = prop_gpts(i,1)   + N(i) * w(integ_point) * dt * damage_el(1,integ_point);            % damage
            prop_gpts(i,2:4) = prop_gpts(i,2:4) + N(i) * w(integ_point) * dt * stress_vec';                         % stresses [sigma_xx,sigma_yy,tau_xy]
            prop_gpts(i,5:7) = prop_gpts(i,5:7) + N(i) * w(integ_point) * dt * [exx,eyy,exy];                       % strains [e_xx, e_yy, gxy]
            prop_gpts(i,8)   = prop_gpts(i,8)   + N(i) * w(integ_point) * dt * e_star;                              % local equivalent strain [e_star]
            prop_gpts(i,9)   = prop_gpts(i,9)   + N(i) * w(integ_point) * dt * nonlocal_strain_el(1,integ_point);   % nonlocal equivalent strain             
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

