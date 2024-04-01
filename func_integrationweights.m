function [w] = func_integrationweights(quad_elem_order)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====================== INTEGRATION POINTS WEIGHTS =======================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Specify just the special case of 2D quadrilateral elements
% -------------------------------------------------------------------------

% 1st order 
if quad_elem_order == 1
    w = [1.,1.,1.,1.]; 

% 2nd order 
elseif quad_elem_order == 2
    w = zeros(1,9);
    w1D = [0.555555555,0.888888888,0.55555555555];
    for j = 1:3
        for i = 1:3
            n    = 3 * (j-1) + i;
            w(n) = w1D(i) * w1D(j);
        end
    end   

else
    disp('Check the quadrilateral element order!')
end


end
