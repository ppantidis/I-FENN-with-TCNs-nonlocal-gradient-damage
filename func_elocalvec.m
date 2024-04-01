function [local_strain_el, nodes_prop_vec_elem] = func_elocalvec(lmncoord,dofs,quad_elem_order,estar_shear_flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======================= ELEMENT LOCAL STRAIN VECTOR =====================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if quad_elem_order == 1
    number_of_nodes = 4;
elseif quad_elem_order == 2
    number_of_nodes = 8;
else
    disp('Check your quad_elem_order variable')
end


% Isolate x-y dofs 
dofs_u = dofs(1:2,:); % 8x1

% Unroll dofs into a vector
dofs_vec = dofs_u(:); % 8x1

% Setting up: 
npoints = func_numberofintegrationpoints(quad_elem_order);   % Number of integration points
xilist  = func_integrationpoints(quad_elem_order);           % Positions of integration points
w       = func_integrationweights(quad_elem_order);          % Weights of integration points

% Initialize matrices and vectors at zero
local_strain_el = zeros(1,npoints);
prop_gpts       = zeros(number_of_nodes,1);
M               = zeros(number_of_nodes,number_of_nodes); 

% Loop over the integration points
for integ_point = 1:npoints

    % ================== SHAPE FUNCTIONS AND DERIVATIVES ==================
    % Compute shape functions && derivatives wrt local coords
    xi      = xilist(:,integ_point);             % Size: 2x1
    N       = func_shapefunctions(xi,quad_elem_order);        % Size: 4x1 
    dNdxi   = func_shapefunctionderivs(xi,quad_elem_order);   % Size: 4x2
    
    % Compute the jacobian matrix && its determinant
    dxdxi = lmncoord * dNdxi; % Size: 2x2 - Jacobian matrix
    dt    = det(dxdxi);

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
    
    % Calculate the equivalent strain e_star and the derivatives {s}
    if estar_shear_flag == 0
        [e_star, s] = func_estar_tension(exx, eyy, gxy);
    elseif estar_shear_flag == 1
        [e_star, s] = func_estar_shear(exx, eyy, gxy);
    else
        disp('Check your estar_shear_flag! - func_elstif_Nonlocgradient')
    end
    local_strain_el(1,integ_point) = e_star;

    % Add NaN condition for [s]
    if isnan(s); s(:) = 0; end

    % -----------------------------------------------------------------
    % ===== ASSEMBLE GAUSS POINT PROPERTIES FOR NODAL PROJECTION ======
    % -----------------------------------------------------------------
    for i = 1:number_of_nodes % number of nodes per element: 4
        prop_gpts(i,1)   = prop_gpts(i,1) + N(i) * w(integ_point) * dt * e_star; % local equivalent strain [e_star]
    end

    % Calculate M matrix for projections
    M = M + (N * N' * w(integ_point) * dt);

end

% ====================== COMPUTE NODAL PROPERTIES =====================
nodes_prop_vec_elem = M \ prop_gpts;

end

