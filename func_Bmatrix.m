function [Bu, Be] = func_Bmatrix(dNdx,quad_elem_order)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ===================== STRAIN-DISPLACEMENT MATRIX (B) ====================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% NON LOCAL GRADIENT IMPLEMENTATION %%%%% 
% Defines the strain-displacement matrix B for each element point, by appropriately expanding the dNdx input argument. 
% Bu (3x8) is the strain-displacement matrix used in the calculation of ux, uy      
% Be (2x4) is the strain-displacement matrix used in the calculation of enonlocal

% -------------------------------------------------------------------------
% Specify just the special case of 2D quadrilateral elements
% -------------------------------------------------------------------------

% Transpose the dNdx matrix (for 2D quadrilateral convert from 4x2 to 2x4)
dNdx_transp = dNdx';

% 1st order 
if quad_elem_order == 1

    Bu = zeros(3,8);
    Bu(1,1) = dNdx_transp(1,1);
    Bu(2,2) = dNdx_transp(2,1);
    Bu(3,1) = dNdx_transp(2,1);
    Bu(3,2) = dNdx_transp(1,1);

    Bu(1,3) = dNdx_transp(1,2);
    Bu(2,4) = dNdx_transp(2,2);
    Bu(3,3) = dNdx_transp(2,2);
    Bu(3,4) = dNdx_transp(1,2);

    Bu(1,5) = dNdx_transp(1,3);
    Bu(2,6) = dNdx_transp(2,3);
    Bu(3,5) = dNdx_transp(2,3);
    Bu(3,6) = dNdx_transp(1,3);

    Bu(1,7) = dNdx_transp(1,4);
    Bu(2,8) = dNdx_transp(2,4);
    Bu(3,7) = dNdx_transp(2,4);
    Bu(3,8) = dNdx_transp(1,4);
    
    Be = dNdx_transp;

% 2nd order 
elseif quad_elem_order == 2

    Bu = zeros(3,16);

    Bu(1,1) = dNdx_transp(1,1);
    Bu(2,2) = dNdx_transp(2,1);
    Bu(3,1) = dNdx_transp(2,1);
    Bu(3,2) = dNdx_transp(1,1);

    Bu(1,3) = dNdx_transp(1,2);
    Bu(2,4) = dNdx_transp(2,2);
    Bu(3,3) = dNdx_transp(2,2);
    Bu(3,4) = dNdx_transp(1,2);

    Bu(1,5) = dNdx_transp(1,3);
    Bu(2,6) = dNdx_transp(2,3);
    Bu(3,5) = dNdx_transp(2,3);
    Bu(3,6) = dNdx_transp(1,3);
    
    Bu(1,7) = dNdx_transp(1,4);
    Bu(2,8) = dNdx_transp(2,4);
    Bu(3,7) = dNdx_transp(2,4);
    Bu(3,8) = dNdx_transp(1,4);

    Bu(1,9)  = dNdx_transp(1,5);
    Bu(2,10) = dNdx_transp(2,5);
    Bu(3,9)  = dNdx_transp(2,5);
    Bu(3,10) = dNdx_transp(1,5);

    Bu(1,11) = dNdx_transp(1,6);
    Bu(2,12) = dNdx_transp(2,6);
    Bu(3,11) = dNdx_transp(2,6);
    Bu(3,12) = dNdx_transp(1,6);

    Bu(1,13) = dNdx_transp(1,7);
    Bu(2,14) = dNdx_transp(2,7);
    Bu(3,13) = dNdx_transp(2,7);
    Bu(3,14) = dNdx_transp(1,7);

    Bu(1,15) = dNdx_transp(1,8);
    Bu(2,16) = dNdx_transp(2,8);
    Bu(3,15) = dNdx_transp(2,8);
    Bu(3,16) = dNdx_transp(1,8);

    Be = dNdx_transp;

else
    disp('Check the quadrilateral element order!')
end

end






