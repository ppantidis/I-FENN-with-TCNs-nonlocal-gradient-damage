function [n] = func_numberofintegrationpoints(quad_elem_order)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======================== NO. INTEGRATION POINTS =========================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Specify just the special case of 2D quadrilateral elements
% -------------------------------------------------------------------------

% 1st order 
if quad_elem_order == 1
    n = 4; 

% 2nd order 
elseif quad_elem_order == 2
    n = 9;
    
else
    disp('Check the quadrilateral element order!')
end

end
