function [xi] = func_integrationpoints(quad_elem_order)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ===================== INTEGRATION POINTS POSITIONS ======================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Specify just the special case of 2D quadrilateral elements
% -------------------------------------------------------------------------

% 1st order 
if quad_elem_order == 1
    xi      = zeros(2,4);
    xi(1,1) = -0.5773502692;
    xi(2,1) = xi(1,1);
    xi(1,2) = -xi(1,1);
    xi(2,2) = xi(1,1);
    xi(1,3) = xi(1,1);
    xi(2,3) = -xi(1,1);
    xi(1,4) = -xi(1,1);
    xi(2,4) = -xi(1,1);

% 2nd order 
elseif quad_elem_order == 2
    xi = zeros(2,9);
    xi(1,1) = -0.7745966692;
    xi(2,1) = xi(1,1);
    xi(1,2) = 0.0;
    xi(2,2) = xi(1,1);
    xi(1,3) = -xi(1,1);
    xi(2,3) = xi(1,1);
    xi(1,4) = xi(1,1);
    xi(2,4) = 0.0;
    xi(1,5) = 0.0;
    xi(2,5) = 0.0;
    xi(1,6) = -xi(1,1);
    xi(2,6) = 0.0;
    xi(1,7) = xi(1,1);
    xi(2,7) = -xi(1,1);
    xi(1,8) = 0.;
    xi(2,8) = -xi(1,1);
    xi(1,9) = -xi(1,1);
    xi(2,9) = -xi(1,1);
    
else
    disp('Check the quadrilateral element order!')
end




end