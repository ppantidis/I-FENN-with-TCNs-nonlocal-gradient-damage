function [d2N_dksi2, d2N_deta2, d2N_dksideta] = func_second_partial_derivs(xi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ==================== PARTIAL DERIVATIVES OF NONLOCAL ====================
% ==================== STRAIN WRT NATURAL COORDINATES =====================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Include global variables
func_include_flags;

d2N_dksi2    = zeros(1,8);
d2N_deta2    = zeros(1,8);
d2N_dksideta = zeros(1,8);

ksi = xi(1); 
eta = xi(2);

% Compute 2nd order derivatives of shape functions wrt ksi
d2N_dksi2(1,1) = (1 - eta)/2;
d2N_dksi2(1,2) = (1 - eta)/2;
d2N_dksi2(1,3) = (1 + eta)/2;
d2N_dksi2(1,4) = (1 + eta)/2;
d2N_dksi2(1,5) = eta - 1;
d2N_dksi2(1,6) = 0;
d2N_dksi2(1,7) = (-eta - 1);
d2N_dksi2(1,8) = 0;


% Compute 2nd order derivatives of shape functions wrt eta
d2N_deta2(1,1) = (1 - ksi)/2;
d2N_deta2(1,2) = (1 + ksi)/2;
d2N_deta2(1,3) = (1 + ksi)/2;
d2N_deta2(1,4) = (1 - ksi)/2;
d2N_deta2(1,5) = 0;
d2N_deta2(1,6) = (-ksi - 1);
d2N_deta2(1,7) = 0;
d2N_deta2(1,8) = (ksi - 1);


% Compute 2nd order derivatives of shape functions wrt ksi and eta
d2N_dksideta(1,1) = (1 - 2*ksi - 2*eta)/4;
d2N_dksideta(1,2) = (2*eta - 2*ksi - 1)/4;
d2N_dksideta(1,3) = (1 + 2*ksi + 2*eta)/4;
d2N_dksideta(1,4) = (2*ksi - 2*eta - 1)/4;
d2N_dksideta(1,5) = ksi;
d2N_dksideta(1,6) = (-eta); 
d2N_dksideta(1,7) = (-ksi);
d2N_dksideta(1,8) = eta;


end








