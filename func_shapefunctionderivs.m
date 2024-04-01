function [dNdxi] = func_shapefunctionderivs(xi,quad_elem_order)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====================== SHAPE FUNCTION DERIVATIVES =======================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculates shape function derivatives for various element types

% -------------------------------------------------------------------------
% Specify just the special case of 2D quadrilateral elements
% -------------------------------------------------------------------------

% 1st order
if quad_elem_order == 1
    dNdxi      = zeros(4,2);
    dNdxi(1,1) = -0.25*(1.-xi(2));
    dNdxi(1,2) = -0.25*(1.-xi(1));
    dNdxi(2,1) = 0.25*(1.-xi(2));
    dNdxi(2,2) = -0.25*(1.+xi(1));
    dNdxi(3,1) = 0.25*(1.+xi(2));
    dNdxi(3,2) = 0.25*(1.+xi(1));
    dNdxi(4,1) = -0.25*(1.+xi(2));
    dNdxi(4,2) = 0.25*(1.-xi(1));

elseif quad_elem_order == 2
    dNdxi       = zeros(8,2);
    dNdxi(1,1)  = 0.25*(1.-xi(2))*(2.*xi(1)+xi(2));
    dNdxi(1,2)  = 0.25*(1.-xi(1))*(xi(1)+2.*xi(2));
    dNdxi(2,1)  = 0.25*(1.-xi(2))*(2.*xi(1)-xi(2));
    dNdxi(2,2)  = 0.25*(1.+xi(1))*(2.*xi(2)-xi(1));
    dNdxi(3,1)  = 0.25*(1.+xi(2))*(2.*xi(1)+xi(2));
    dNdxi(3,2)  = 0.25*(1.+xi(1))*(2.*xi(2)+xi(1));
    dNdxi(4,1)  = 0.25*(1.+xi(2))*(2.*xi(1)-xi(2));
    dNdxi(4,2)  = 0.25*(1.-xi(1))*(2.*xi(2)-xi(1));
    dNdxi(5,1)  = -xi(1)*(1.-xi(2));
    dNdxi(5,2)  = -0.5*(1.-xi(1)*xi(1));
    dNdxi(6,1)  = 0.5*(1.-xi(2)*xi(2));
    dNdxi(6,2)  = -(1.+xi(1))*xi(2);
    dNdxi(7,1)  = -xi(1)*(1.+xi(2));
    dNdxi(7,2)  = 0.5*(1.-xi(1)*xi(1));
    dNdxi(8,1)  = -0.5*(1.-xi(2)*xi(2));
    dNdxi(8,2)  = -(1.-xi(1))*xi(2);
    
else
    disp('Check the quadrilateral element order!')
end

end