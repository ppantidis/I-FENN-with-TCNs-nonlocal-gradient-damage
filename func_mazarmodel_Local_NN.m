function [omega,domega_dkappa] = func_mazarmodel_Local_NN(enonlocal_strain,alpha_val,beta_val,e_delta,dmax)
%================= MAZAR'S DAMAGE MODEL ================================

% Computes the damage variable delta based on:
% e_star      : equivalent strain
% e_delta     : threshold strain
% alpha_val   : model hyperparameter
% beta_val    : model hyperparameter

if enonlocal_strain < e_delta
    omega = 0;
    domega_dkappa = 0;
else
    omega = dmax * (1 - ((e_delta * (1 - alpha_val) / enonlocal_strain) + (alpha_val / exp(beta_val * (enonlocal_strain - e_delta)))));
    domega_dkappa = dmax * (alpha_val * beta_val * exp(beta_val * (e_delta - enonlocal_strain)) + e_delta * (1 - alpha_val) / (enonlocal_strain^2));
end

end