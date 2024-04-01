function [Partial_nonlocal_natural] = func_partial_nonlocal_natural(dofs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ==================== PARTIAL DERIVATIVES OF NONLOCAL ====================
% ==================== STRAIN WRT NATURAL COORDINATES =====================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Include global variables
func_include_flags;

% Setting up: 
npoints = func_numberofintegrationpoints(quad_elem_order);   % Number of integration points   
xilist  = func_integrationpoints(quad_elem_order);           % Positions of integration points

GP_counter = 0;
Partial_nonlocal_natural = zeros(3,1,npoints*nelem);
 
for lmn = 1:nelem
    
    lmncoord = zeros(ncoord,maxnodes);
    % Extract nodal coordinates and dofs for the current element
    for a = 1:nelnodes(lmn)
        for i = 1:ncoord
            lmncoord(i,a) = coords(i,connect(a,lmn));
        end        
        for i = 1:ndof
            lmndof(i,a) = dofs(ndof*(connect(a,lmn)-1)+i);
        end
    end
    
    xcoordinates = lmncoord(1,:)';   % Size: 8x1
    ycoordinates = lmncoord(2,:)';   % Size: 8x1 
    dofs_e       = lmndof(3,:);        % Size: 1x4
    dofs_e_vec   = dofs_e(:);

    for integ_point = 1:npoints

        GP_counter  = GP_counter + 1;                               % Increase counter
        xi          = xilist(:,integ_point);                        % Size: 2x1
        dNdxi       = func_shapefunctionderivs(xi,quad_elem_order); % Size: 4x2
        dxdxi       = lmncoord * dNdxi;                             % Size: 2x2 - Jacobian matrix
        dxidx       = inv(dxdxi);                                   % Size: 2x2 - inverse of Jacobian matrix

        % Note for notation
        % dxidx = [a b; c d], where: 
        a = dxidx(1,1); % a = partial(xi)/partial(x)
        b = dxidx(1,2); % b = partial(eta)/partial(x)
        c = dxidx(2,1); % c = partial(xi)/partial(y)
        d = dxidx(2,2); % d = partial(eta)/partial(y)

        % Compute 2nd order partial derivatives of shape functions wrt
        % natural (d2N/dksi2, d2N/deta2, d2N/dksideta) - each vector is 1x8 
        [d2N_dksi2, d2N_deta2, d2N_dksideta] = func_second_partial_derivs(xi);
      
        % -----------------------------------------------------------------
        % -------------------------- GREEN TERMS --------------------------
        % -----------------------------------------------------------------

        % Compute 2nd order partial derivatives of physical wrt natural 
        % (d2x/dksi2, d2x/deta2, d2y/dksi2, d2y/deta2, d2x/dksideta,
        % d2y/dksideta) - each value is a scalar
        
        d2x_dksi2       = d2N_dksi2     * xcoordinates;
        d2x_deta2       = d2N_deta2     * xcoordinates;
        d2x_dksideta    = d2N_dksideta  * xcoordinates;
        d2y_dksi2       = d2N_dksi2     * ycoordinates;        
        d2y_deta2       = d2N_deta2     * ycoordinates;
        d2y_dksideta    = d2N_dksideta  * ycoordinates;        

        % -----------------------------------------------------------------
        % -------------------------- RED TERMS --------------------------
        % -----------------------------------------------------------------

        % Compute 2nd order partial derivatives of nonlocal strain wrt natural 
        % (d2e/dksi2, d2e/deta2, d2e/dksideta) - each value is a scalar
        
        d2e_dksi2    = d2N_dksi2    * dofs_e_vec;
        d2e_deta2    = d2N_deta2    * dofs_e_vec;
        d2e_dksideta = d2N_dksideta * dofs_e_vec;        


        % -----------------------------------------------------------------
        % -------------------------- BLUE TERMS ---------------------------
        % -----------------------------------------------------------------

        % Compute 1st order partial derivatives of nonlocal strain wrt
        % physical (de/dx, de/dy) - each value is a scalar
        
        de_dx = dNdxi(:,1)' * dofs_e_vec * a + dNdxi(:,2)' * dofs_e_vec * b;
        de_dy = dNdxi(:,1)' * dofs_e_vec * c + dNdxi(:,2)' * dofs_e_vec * d;

        
        % -----------------------------------------------------------------
        % ------------------- PARTIAL NONLOCAL NATURAL --------------------
        % -----------------------------------------------------------------

        Partial_nonlocal_natural(:,1,GP_counter) = [d2e_dksi2    - de_dx * d2x_dksi2    - de_dy * d2y_dksi2; ...
                                                    d2e_deta2    - de_dx * d2x_deta2    - de_dy * d2y_deta2; ...
                                                    d2e_dksideta - de_dx * d2x_dksideta - de_dy * d2y_dksideta];

    end

end


end