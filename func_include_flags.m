%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======================= DEFINE GLOBAL VARIABLES =========================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FEM MODEL PARAMETERS 
global nprops materialprops ncoord ndof nnodes coords nelem maxnodes connect nelnodes elident_vec nfix fixnodes ndof1 ndof2 quad_elem_order
ndof1 = 1;
ndof2 = 2;

% NONLOCAL GRADIENT PARAMETER 
% g = lc^2/2, lc = characteristic length
global g

% MAZAR'S DAMAGE MODEL PARAMETERS
global alpha_val beta_val e_delta dmax estar_shear_flag

% ADAPTIVE LOAD AND PLOTTING PARAMETERS
global min_iter max_iter max_accept_iter dlfactor_incr_threshold increment_plot_threshold loadfactor_plot_threshold

% SOLVER SCHEME
% SolverID: 1 - Local, 2 - Nonlocal Gradient, 3 - Nonlocal Integral
% TangentID: 1 - Analytical, 2 - Numerical
global SolverID TangentID
