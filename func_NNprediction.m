%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ============== NEURAL NETWORK PREDICTION OF NONLOCAL STRAIN =============
% ================== AND ITS DERIVATIVE WRT LOCAL STRAIN ==================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% The following code goes inside the MAIN_Script

fid = py.open("Trained_model_xNOTnormalized_2hidlay_2units_10ep.pickle",'rb');
NNweights = py.pickle.load(fid);
NNweights = cell(NNweights)

xcoord = 1.3333;
ycoord = 0.42265;
g_val  = 8;
elocal = 0.000000285;
[enonlocal, denonlocal_dlocal] = func_NNprediction(xcoord,ycoord,g_val,elocal,NNweights);


% -------------------------------------------------------------------------
% The following code stays inside the func_NNprediction

function [enonlocal, denonlocal_dlocal] = func_NNprediction(xcoord,ycoord,g,elocal,NNweights)

X = [xcoord ycoord elocal g]'

for i = 1:length(NNweights{1})-1
    W = double(NNweights{1}{i})'
    b = double(NNweights{2}{i})'
    X = tanh(W*X + b)
end
W = double(NNweights{1}{length(NNweights{1})})'
b = double(NNweights{2}{length(NNweights{1})})'
enonlocal = W*X + b


denonlocal_dlocal = 0




end










