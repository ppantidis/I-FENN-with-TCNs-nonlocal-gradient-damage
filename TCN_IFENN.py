### Import packages ###
import torch #, torchvision.models
import torch.nn as nn
import torch.nn.functional as F
from tcn import TemporalConvNet as TCN
import numpy as np
import scipy.io
import os 
import time
import pickle
import numpy.random as npr
import random
import warnings
import gc

gc.collect()
warnings.filterwarnings('ignore')
torch.cuda.empty_cache()
PYTORCH_CUDA_ALLOC_CONF = "max_split_size_mb:<128>"

###########################################################################################################################################
# --------------------------------------------------------------------------------
Model_name      = "SNS_Coarse_D10_P00_dec10000_Ad2000_LBFGS5000"
fileDir         = Model_name + ".pickle"

# --------------------------------------------------------------------------------
start_inc_train   = 1
step_inc_train    = 1
end_inc_train     = 160
Nincr_train       = int((end_inc_train - start_inc_train) / step_inc_train + 1)

start_inc_pred  = 1
step_inc_pred   = 1
end_inc_pred    = 160
Nincr_pred      = int((end_inc_pred - start_inc_pred) / step_inc_pred + 1)

# Assign the value to the g variable (g = lc^2/2)
g = 4.5

###########################################################################################################################################
# --------------------------------------------------------------------------------
Function_flag       = 0                 # Options: {1 (Train mode), 0 (Predict mode)}
Cost_flag           = "data"            # Options: {"data", "physics", "data_physics"}
Scaling_flag        = "decimal"         # Options: {"decimal", "amsf", "minmax"}
Optimization_flag   = "SF"            # Options: {"SF" (shape functions), "AD" (automatic differentiation)}

# --------------------------------------------------------------------------------
decimal_scaling       = 10000                                                                           # Decimal scaling factor, options: {1 (no scaling), 100, 1000, 10000}
amsf_scaling          = np.array(list(scipy.io.loadmat('amsf.mat').items()), dtype = object)[3,1]   # Load amsf scaling factors
minmax_upper_bound    = 10                                                                             # Upper bound in min-max range, options: {[0 1], [0 10], [0 100]}

###########################################################################################################################################
# --------------------------------------------------------------------------------
Dimensions_input    = 4 # input (X): x-coord, y-coord, elocal, loadfactor 
Dimensions_output   = 1 # output (Y): enonlocal
num_tcn_layers      = 1
num_filters         = 6
num_dilations       = 3
num_channels1       = [num_filters] * num_dilations 
num_channels2       = [num_filters] * num_dilations
kernel_size         = 12
dropout             = 0.00
activ_func          = "tanh"
mesh_file           = "SNS_Coarse"

# --------------------------------------------------------------------------------
Adam_epochs  = 2000                  # Number of Adam epochs
Adam_lr      = 0.001                 # Adam learning rate
LBFGS_epochs = 5000                # Number of L-BFGS epochs

# --------------------------------------------------------------------------------
Adam_print_epochs_inc   = 100         # Print cost value every no. epochs of Adam
LBFGS_print_epochs_inc  = 100         # Print cost value every no. epochs of L-BFGS
Adam_save_epochs_inc    = 50000         # Save trained model every no. epochs of Adam
LBFGS_save_epochs_inc   = 50000         # Save trained model every no. epochs of L-BFGS

###########################################################################################################################################
# Fix seed 
seeding_number = 1; npr.seed(seeding_number); torch.manual_seed(seeding_number); random.seed(seeding_number)

# Define printing precision
torch.set_printoptions(precision = 5)

###########################################################################################################################################
# --------------------------------------------------------------------------------
Enonlocal_term_norm_vec = []    # Vector to store the L2 norm of the enonlocal term in the residual PDE  
Laplacian_term_norm_vec = []    # Vector to store the L2 norm of the laplacian term in the residual PDE  
Elocal_term_norm_vec = []       # Vector to store the L2 norm of the elocal term in the residual PDE 
loss_vec = []                   # Vector to store the total loss  
loss_PDE_GP_vec = []            # Vector to store the PDE loss   
loss_PDEBCs_nodeslrb_vec = []   # Vector to store the PDEBCs loss at left-right boundary  
loss_PDEBCs_nodesbtb_vec = []   # Vector to store the PDEBCs loss at bottom-top boundary
loss_Data_GP_vec = []           # Vector to store the data loss  
Iter_LBFGS = [0]                # Scalar to store the number of LBFGS iterations 
Training_times_vec = []         # Vector to store the training times (Adam, L-BFGS, total)

# --------------------------------------------------------------------------------
# Note: Including the data and/or the physics is also accounted for in the cost function definition. 
if Cost_flag == "data":
    w_PDE_GP, w_PDEBCs_nodeslrb, w_PDEBCs_nodesbtb = 0, 0, 0    # Weight of the PDE loss term (at Gauss Points) and PDEBCs loss term (at nodes in the left-right boundary & bottom-top boundary) 
    w_Data_GP = 1                                               # Weight of the labeled data loss term (at Gauss Points)
elif Cost_flag == "physics":
    w_PDE_GP, w_PDEBCs_nodeslrb, w_PDEBCs_nodesbtb = 1, 1, 1 
    w_Data_GP = 0
elif Cost_flag == "data_physics":
    w_PDE_GP, w_PDEBCs_nodeslrb, w_PDEBCs_nodesbtb = 0.5, 0.5, 0.5
    w_Data_GP = 0.5
else:
    print("Check your Cost_flag variable!") 

###########################################################################################################################################
if Function_flag == 1:

    # ----------------------------------------------------------------------------
    ################################# Train mode #################################
    # ----------------------------------------------------------------------------

    device = torch.device('cuda')  # Options: {'cuda', 'cpu'}
    Nincr = Nincr_train
 
    # Load training data 
    All_TrainingData_SNS_Coarse_GP_np = np.array(list(scipy.io.loadmat('TrainingData_SNS_Coarse_GP_time_features_lf_conv1e4.mat').items()), dtype = object)[3,1]             

    # Slice training data
    All_TrainingData_SNS_Coarse_GP_np = All_TrainingData_SNS_Coarse_GP_np[:,start_inc_train-1:end_inc_train:step_inc_train,:]        # GPs   x increments x features

    # Scale training data
    if Scaling_flag == "decimal":
        All_TrainingData_SNS_Coarse_GP_np[:,:,2] = np.multiply(All_TrainingData_SNS_Coarse_GP_np[:,:,2], decimal_scaling)        
        All_TrainingData_SNS_Coarse_GP_np[:,:,3] = np.multiply(All_TrainingData_SNS_Coarse_GP_np[:,:,3], decimal_scaling)  
        
    elif Scaling_flag == "amsf":

        amsf_scaling = amsf_scaling[0,start_inc_train-1:end_inc_train:step_inc_train]

        All_TrainingData_SNS_Coarse_GP_np[:,:,2] = np.multiply(All_TrainingData_SNS_Coarse_GP_np[:,:,2], amsf_scaling)        
        All_TrainingData_SNS_Coarse_GP_np[:,:,3] = np.multiply(All_TrainingData_SNS_Coarse_GP_np[:,:,3], amsf_scaling)  

    elif Scaling_flag == "minmax":

        elocal_min  = np.min(All_TrainingData_SNS_Coarse_GP_np[:,:,2], axis = 0) # Calculate elocal min at each time increment
        elocal_max  = np.max(All_TrainingData_SNS_Coarse_GP_np[:,:,2], axis = 0) # Calculate elocal max at each time increment        
        a           = minmax_upper_bound / (elocal_max - elocal_min)
        b           = - (minmax_upper_bound * elocal_min) / (elocal_max - elocal_min)        

        All_TrainingData_SNS_Coarse_GP_np[:,:,2] = a * All_TrainingData_SNS_Coarse_GP_np[:,:,2] + b
        All_TrainingData_SNS_Coarse_GP_np[:,:,3] = a * All_TrainingData_SNS_Coarse_GP_np[:,:,3] + b
        
    else:
        print("Check your Scaling_flag variable!") 

    # Check for NaN or inf values in the datasets
    if np.isnan(All_TrainingData_SNS_Coarse_GP_np).any() or np.isinf(All_TrainingData_SNS_Coarse_GP_np).any():
        print("All_TrainingData_SNS_Coarse_GP_np contains NaN or infinite values!")

    # Placeholder for predict data
    PredictData_SNS_Coarse_GP_np = []

elif Function_flag == 0:

    # ----------------------------------------------------------------------------
    ################################ Predict mode ################################
    # ----------------------------------------------------------------------------

    device = torch.device('cuda')  # Options: {'cuda', 'cpu'}
    Nincr = Nincr_pred

    # Note: "increment" and "iteration" are provided during I-FENN execution
    TrainingData_placeholder_GP_np = np.array(list(scipy.io.loadmat("TrainingData_SNS_Coarse_GP_time_features_lf_conv1e4.mat").items()), dtype = object)[3,1]

    # Load prediction data 
    # PredictData_SNS_Coarse_GP_np = np.array(list(scipy.io.loadmat("TrainingData_SNS_Coarse_GP_time_features_lf.mat").items()), dtype = object)[3,1]
    PredictData_SNS_Coarse_GP_np = np.array(list(scipy.io.loadmat("PredData_IFENN_GP_time_features_lf_" + str(increment) + "_" + str(iteration) + ".mat").items()), dtype = object)[3,1]

    # Scale prediction data
    if Scaling_flag == "decimal":
        PredictData_SNS_Coarse_GP_np[:,:,2] = np.multiply(PredictData_SNS_Coarse_GP_np[:,:,2], decimal_scaling)        
     
    elif Scaling_flag == "amsf":
        amsf_scaling_predict = amsf_scaling[0,int(increment)-1]
        PredictData_SNS_Coarse_GP_np[:,:,2] = np.multiply(PredictData_SNS_Coarse_GP_np[:,:,2], amsf_scaling_predict)    
     
    elif Scaling_flag == "minmax":

        elocal_min  = np.min(TrainingData_placeholder_GP_np[:,:,2], axis = 0) # Calculate elocal min at each time increment
        elocal_max  = np.max(TrainingData_placeholder_GP_np[:,:,2], axis = 0) # Calculate elocal max at each time increment
        a           = minmax_upper_bound / (elocal_max - elocal_min)
        b           = - (minmax_upper_bound * elocal_min) / (elocal_max - elocal_min)

        PredictData_SNS_Coarse_GP_np[:,:,2] = (a * PredictData_SNS_Coarse_GP_np[:,:,2] + b)
     
    else:
        print("Check your Scaling_flag variable!")

    PredictData_SNS_Coarse_GP_np = np.nan_to_num(PredictData_SNS_Coarse_GP_np)    
    
    # Placeholder for train data
    All_TrainingData_SNS_Coarse_GP_np = [], [], []

else:
    print("Check your Function_flag variable!") 


###########################################################################################################################################
# Build TCN class 
class Seq2Seq(nn.Module):
    def __init__(self, Dimensions_input, Dimensions_output, PredictData_SNS_Coarse_GP_np, All_TrainingData_SNS_Coarse_GP_np):
        super(Seq2Seq, self).__init__()

        # ----------------------------------------------------------------------------------------------------------------        
        # Convert input data to tensors         
        if Function_flag == 1:
            self.X_TrainingData_SNS_Coarse_GP = torch.tensor(All_TrainingData_SNS_Coarse_GP_np[:,:,[0, 1, 2, 4]], requires_grad = True).float().to(device)
            self.Y_TrainingData_SNS_Coarse_GP = torch.tensor(All_TrainingData_SNS_Coarse_GP_np[:,:,[3]], requires_grad = False).float().to(device)
        elif Function_flag == 0:
            self.X_PredictData_SNS_Coarse_GP = torch.tensor(PredictData_SNS_Coarse_GP_np[:,:,[0, 1, 2, 4]], requires_grad = True).float().to(device)
        else:
            print("Check your Function_flag variable!") 

        # ----------------------------------------------------------------------------------------------------------------
        # Define TCN layers
        self.tcn1   = TCN(Dimensions_input, num_channels1, activ_func, kernel_size, dropout).to(device)
        self.tcn2   = TCN(num_channels1[-1], num_channels2, activ_func, kernel_size, dropout).to(device)
        self.linear1 = nn.Linear(num_filters, Dimensions_output).to(device)

        # ----------------------------------------------------------------------------------------------------------------
        # Define optimizers
        self.optimizer_LBFGS = torch.optim.LBFGS(
            self.parameters(), 
            lr                  = 1.6, 
            max_iter            = LBFGS_epochs, 
            max_eval            = LBFGS_epochs * 1.25, 
            history_size        = 100,
            tolerance_grad      = 1e-9, 
            tolerance_change    = 1.0 * np.finfo(float).eps,
            line_search_fn      = "strong_wolfe")
        
        self.optimizer_Adam = torch.optim.Adam(self.parameters(), lr = Adam_lr)
        self.iter = 0

        # Initialize NNs
        if Function_flag == 1:
            self.parameters = self.init_weights()

    # ----------------------------------------------------------------------------------------------------------------
    # Make enonlocal prediction
    def forward(self, x): 

        # 1 TCN layer
        if num_tcn_layers == 1:
            y1  = self.tcn1(x.transpose(1,2))
            out = self.linear1(y1.transpose(1,2))

        # 2 TCN layers
        elif num_tcn_layers == 2:
            y1  = self.tcn1(x.transpose(1,2))
            y2  = self.tcn2(y1)
            out = self.linear1(y2.transpose(1,2))

        else:
            print("Check your num_to_layers flag!")

        return out

    # ----------------------------------------------------------------------------------------------------------------
    # Calculate 1st order partial derivative of enonlocal (combo of "enonlocal_elocal" and "enonlocal_gradient")
    def first_partial_enonlocal(self, enonlocal, X):

        # Compute 1st order derivatives wrt the input features: x-coord, y-coord, elocal, lf
        denonlocal_dinput = torch.autograd.grad(enonlocal, X, grad_outputs=torch.ones_like(enonlocal), create_graph=True, retain_graph=True)[0].to(device)

        # Initialize empty tensors
        enonlocal_grad_xy = torch.empty((len(X), Nincr, 2)).to(device)
        enonlocal_elocal  = torch.empty((len(X), Nincr, 1)).to(device)

        enonlocal_grad_xy[:,:,0] = denonlocal_dinput[:, :, 0].to(device) # w.r.t. x_coord
        enonlocal_grad_xy[:,:,1] = denonlocal_dinput[:, :, 1].to(device) # w.r.t. y_coord
        enonlocal_elocal[:,:,0]  = denonlocal_dinput[:, :, 2].to(device) # w.r.t. elocal
        
        return enonlocal_grad_xy, enonlocal_elocal

    # ----------------------------------------------------------------------------------------------------------------
    # Calculate 1st order partial derivative of enonlocal wrt elocal
    def enonlocal_elocal(self, enonlocal, X):

        # input features: x-coord, y-coord, elocal, lf
        denonlocal_dinput = torch.autograd.grad(enonlocal, X, grad_outputs=torch.ones_like(enonlocal), create_graph=True, retain_graph=True)[0].to(device)
        enonlocal_elocal = torch.empty((len(X), Nincr, 1)).to(device)
        enonlocal_elocal[:,:,0] = denonlocal_dinput[:, :, 2].to(device) # w.r.t. elocal
        
        return enonlocal_elocal
    
    # ----------------------------------------------------------------------------------------------------------------
    # Calculate 1st order partial derivative of enonlocal wrt x-coord and y-coord
    def enonlocal_gradient(self, enonlocal, X):

        # input features: x-coord, y-coord, elocal, lf
        denonlocal_dinput = torch.autograd.grad(enonlocal, X, grad_outputs=torch.ones_like(enonlocal), create_graph=True, retain_graph=True)[0].to(device)    
        enonlocal_grad_xy = torch.empty((len(X), Nincr, 2)).to(device)
        enonlocal_grad_xy[:,:,0] = denonlocal_dinput[:, :, 0].to(device) # w.r.t. x_coord
        enonlocal_grad_xy[:,:,1] = denonlocal_dinput[:, :, 1].to(device) # w.r.t. y_coord

        return enonlocal_grad_xy

    # ----------------------------------------------------------------------------------------------------------------
    # Calculate 2nd order partial derivative of enonlocal wrt x-coord and y-coord
    def enonlocal_laplacian(self, enonlocal_grad_xy, X):

        # input features: x-coord, y-coord, elocal, lf
        dq1dxy = torch.autograd.grad(enonlocal_grad_xy[:,:,0], X, grad_outputs=torch.ones_like(enonlocal_grad_xy[:,:,0]), create_graph=True, retain_graph=True)[0]
        dq2dxy = torch.autograd.grad(enonlocal_grad_xy[:,:,1], X, grad_outputs=torch.ones_like(enonlocal_grad_xy[:,:,1]), create_graph=True, retain_graph=True)[0]
        dq1dx  = dq1dxy[:, :, 0].unsqueeze(2) # w.r.t. x_coord
        dq2dy  = dq2dxy[:, :, 1].unsqueeze(2) # w.r.t. y_coord
        enonlocal_laplacian = dq1dx + dq2dy
        return enonlocal_laplacian    
    
    # ----------------------------------------------------------------------------------------------------------------
    # Define the LBFGS loss function
    def LBFGS_loss_function(self):
        
        # Make prediction
        enonlocal_pred_GP       = self.forward(self.X_TrainingData_SNS_Coarse_GP)

        # ---------------------------------------------------------------------------------------------------------------
        if Cost_flag == "data" or Cost_flag == "data_physics":
            Loss_Data_GP = torch.squeeze(torch.linalg.vector_norm(enonlocal_pred_GP - self.Y_TrainingData_SNS_Coarse_GP, ord = 2, dim = 0)).to(device)
            
            # Placheholder for the PDE loss terms 
            if Cost_flag == "data":
                Loss_PDE_GP, Loss_PDEBCs_nodeslrb, Loss_PDEBCs_nodesbtb = torch.zeros(Nincr).to(device), torch.zeros(Nincr).to(device), torch.zeros(Nincr).to(device)         

        # ---------------------------------------------------------------------------------------------------------------
        if Cost_flag == "physics" or Cost_flag == "data_physics":
            
            # Computations at Gauss Points
            if Optimization_flag == "AD":
                # Calculate 1st order partial derivatives
                enonlocal_grad_xy_GP = self.enonlocal_gradient(enonlocal_pred_GP, self.X_TrainingData_SNS_Coarse_GP)
                # Calculate 2nd order partial derivatives
                laplacian_enonlocal_GP = self.enonlocal_laplacian(enonlocal_grad_xy_GP, self.X_TrainingData_SNS_Coarse_GP)

                # Define the scaling of the laplacian term
                a_scaling = 1                
            
            elif Optimization_flag == "SF":

                # Define the scaling of the laplacian term
                if Scaling_flag == "decimal":
                    a_scaling = torch.unsqueeze(torch.tensor(np.int64(decimal_scaling).astype(np.float32)), dim=-1).to(device)                           
                elif Scaling_flag == "amsf":
                    a_scaling = torch.unsqueeze(torch.tensor(amsf_scaling.astype(np.float32)), dim=-1).to(device)
                elif Scaling_flag == "minmax":
                    a_scaling = torch.unsqueeze(torch.tensor(a.astype(np.float32)), dim=-1).to(device)
                
            else:
                print("Check your Optimization flag variable!")     

            # Calculate the PDE residuals
            Residual_PDE_GP             = enonlocal_pred_GP - a_scaling * g * laplacian_enonlocal_GP - self.X_TrainingData_SNS_Coarse_GP[:,:,2].unsqueeze(2) 

            # Compute the 2nd norm of the PDE residual terms
            Enonlocal_term_norm     = torch.linalg.vector_norm(enonlocal_pred_GP, ord = 2)
            Laplacian_term_norm     = torch.linalg.vector_norm(a_scaling * g * laplacian_enonlocal_GP, ord = 2)
            Elocal_term_norm        = torch.linalg.vector_norm(self.X_TrainingData_SNS_Coarse_GP[:,:,2].unsqueeze(2), ord = 2)        

            # Compute loss components 
            Loss_PDE_GP             = torch.squeeze(torch.linalg.vector_norm(Residual_PDE_GP, ord = 2, dim = 0)).to(device)
            
            # Store cost values at the end of L-BFGS
            Enonlocal_term_norm_vec.append(Enonlocal_term_norm.detach().cpu().numpy())
            Laplacian_term_norm_vec.append(Laplacian_term_norm.detach().cpu().numpy())
            Elocal_term_norm_vec.append(Elocal_term_norm.detach().cpu().numpy())  

        # Placheholder for the data loss term 
        if Cost_flag == "physics":
            Loss_Data_GP = torch.zeros(Nincr).to(device)

        # ---------------------------------------------------------------------------------------------------------------
        # Compute total loss (contribution of data, boundary nodes, regularization...)
        loss = torch.sum(w_Data_GP * Loss_Data_GP + w_PDE_GP * Loss_PDE_GP + w_PDEBCs_nodeslrb * Loss_PDEBCs_nodeslrb + w_PDEBCs_nodesbtb * Loss_PDEBCs_nodesbtb)
                
        loss_vec.append(loss.detach().cpu().numpy())                                                    
        loss_PDE_GP_vec.append(Loss_PDE_GP.detach().cpu().numpy())
        loss_PDEBCs_nodeslrb_vec.append(Loss_PDEBCs_nodeslrb.detach().cpu().numpy()) 
        loss_PDEBCs_nodesbtb_vec.append(Loss_PDEBCs_nodesbtb.detach().cpu().numpy()) 
        loss_Data_GP_vec.append(Loss_Data_GP.detach().cpu().numpy())

        # Backward and optimize
        self.optimizer_LBFGS.zero_grad()
        loss.backward()   

        self.iter += 1
        if self.iter % LBFGS_print_epochs_inc == 0:
            print('It: %d, Loss: %e ' % (self.iter, loss.item()))  # Print cost values at L-BFGS 

        if self.iter % LBFGS_save_epochs_inc == 0:
            predictions         = self.forward(self.X_TrainingData_SNS_Coarse_GP)
            scipy.io.savemat('Predictions_LBFGS_' + str(self.iter) + '_' + Model_name + '.mat', 
                            {'predictions': predictions.detach().cpu().numpy()})

        Iter_LBFGS[0] = self.iter  
        return loss

    # ----------------------------------------------------------------------------------------------------------------
    def train(self, nIter):
        
        # Adam initiates
        print("--------------------")
        print("Adam training begins")
        for epoch in range(nIter):

            # Make prediction
            enonlocal_pred_GP       = self.forward(self.X_TrainingData_SNS_Coarse_GP)           

            # ---------------------------------------------------------------------------------------------------------------
            if Cost_flag == "data" or Cost_flag == "data_physics":
                Loss_Data_GP = torch.squeeze(torch.linalg.vector_norm(enonlocal_pred_GP - self.Y_TrainingData_SNS_Coarse_GP, ord = 2, dim = 0)).to(device)       
                                
                # Placheholder for the PDE loss terms 
                if Cost_flag == "data":
                    Loss_PDE_GP, Loss_PDEBCs_nodeslrb, Loss_PDEBCs_nodesbtb = torch.zeros(Nincr).to(device), torch.zeros(Nincr).to(device), torch.zeros(Nincr).to(device)  
                    
            # ---------------------------------------------------------------------------------------------------------------
            if Cost_flag == "physics" or Cost_flag == "data_physics":

                # Computations at Gauss Points
                if Optimization_flag == "AD":
                    # Calculate 1st order partial derivatives
                    enonlocal_grad_xy_GP = self.enonlocal_gradient(enonlocal_pred_GP, self.X_TrainingData_SNS_Coarse_GP)

                    # Calculate 2nd order partial derivatives
                    laplacian_enonlocal_GP = self.enonlocal_laplacian(enonlocal_grad_xy_GP, self.X_TrainingData_SNS_Coarse_GP)

                    # Define the scaling of the laplacian term
                    a_scaling = 1

                elif Optimization_flag == "SF":

                    # Define the scaling of the laplacian term
                    if Scaling_flag == "decimal":
                        a_scaling = torch.unsqueeze(torch.tensor(np.int64(decimal_scaling).astype(np.float32)), dim=-1).to(device)                    
                    elif Scaling_flag == "amsf":
                        a_scaling = torch.unsqueeze(torch.tensor(amsf_scaling.astype(np.float32)), dim=-1).to(device)
                    elif Scaling_flag == "minmax":
                        a_scaling = torch.unsqueeze(torch.tensor(a.astype(np.float32)), dim=-1).to(device)
         
                else:
                    print("Check your Optimization flag variable!")

                # Compute the PDE residual
                Residual_PDE_GP             = enonlocal_pred_GP - a_scaling * g * laplacian_enonlocal_GP - self.X_TrainingData_SNS_Coarse_GP[:,:,2].unsqueeze(2)

                # Compute the 2nd norm of the PDE residual terms
                Enonlocal_term_norm     = torch.linalg.vector_norm(enonlocal_pred_GP, ord = 2)
                Laplacian_term_norm     = torch.linalg.vector_norm(a_scaling * g * laplacian_enonlocal_GP, ord = 2)
                Elocal_term_norm        = torch.linalg.vector_norm(self.X_TrainingData_SNS_Coarse_GP[:,:,2].unsqueeze(2), ord = 2)        

                # Compute loss components 
                Loss_PDE_GP             = torch.squeeze(torch.linalg.vector_norm(Residual_PDE_GP, ord = 2, dim = 0)).to(device)
                
                # Store cost values at the end of Adam
                Enonlocal_term_norm_vec.append(Enonlocal_term_norm.detach().cpu().numpy())
                Laplacian_term_norm_vec.append(Laplacian_term_norm.detach().cpu().numpy())            
                Elocal_term_norm_vec.append(Elocal_term_norm.detach().cpu().numpy())

                # Placheholder for the data loss term 
                if Cost_flag == "physics":
                    Loss_Data_GP = torch.zeros(Nincr).to(device)

            # ---------------------------------------------------------------------------------------------------------------
            # Compute total loss (contribution of data, boundary nodes, regularization...)
            loss = torch.sum(w_Data_GP * Loss_Data_GP + w_PDE_GP * Loss_PDE_GP + w_PDEBCs_nodeslrb * Loss_PDEBCs_nodeslrb + w_PDEBCs_nodesbtb * Loss_PDEBCs_nodesbtb)

            loss_vec.append(loss.detach().cpu().numpy())                                                    
            loss_PDE_GP_vec.append(Loss_PDE_GP.detach().cpu().numpy())
            loss_PDEBCs_nodeslrb_vec.append(Loss_PDEBCs_nodeslrb.detach().cpu().numpy()) 
            loss_PDEBCs_nodesbtb_vec.append(Loss_PDEBCs_nodesbtb.detach().cpu().numpy()) 
            loss_Data_GP_vec.append(Loss_Data_GP.detach().cpu().numpy())

            # Backward and optimize
            self.optimizer_Adam.zero_grad()
            loss.backward()
            self.optimizer_Adam.step()
            
            self.iter += 1
            if self.iter % Adam_print_epochs_inc == 0:
                print('It: %d, Loss: %e' %  (self.iter, loss.sum().item()))    # Print cost values at Adam 
                
            if self.iter % Adam_save_epochs_inc == 0:
                predictions         = self.forward(self.X_TrainingData_SNS_Coarse_GP)
                scipy.io.savemat('Predictions_Adam_' + str(self.iter) + '_' + Model_name + '.mat', 
                                {'predictions': predictions.detach().cpu().numpy()})

        Adam_training_time = time.time() - Adam_start_time
        Training_times_vec.append(Adam_training_time)

        # LBFGS initiates
        print("----------------------")
        print("L-BFGS training begins")
        LBFGS_start_time = time.time()
        self.iter = 0
        self.optimizer_LBFGS.step(self.LBFGS_loss_function)

        LBFGS_training_time = time.time() - LBFGS_start_time
        Training_times_vec.append(LBFGS_training_time)

        path = ('./' + Model_name + '.pickle')
        torch.save(model, path)


    # ----------------------------------------------------------------------------------------------------------------
    # Save the trained network
    def save_NN(self, fileDir):
        with open(fileDir, 'wb') as f:
            torch.save([self.parameters.state_dict()], f)
            print("Saved NN parameters successfully!")
        print("----------------------")


    # ---------------------------------------------------------------------------------------------------------------------------------------
    # Load an existing network
    def load_NN(self, fileDir):
        self.parameters = []
        with open(fileDir, 'rb') as f:
            self.parameters = pickle.load(f)
        return self.parameters


    # ----------------------------------------------------------------------------------------------------------------
    # Initialize the network weights
    def init_weights(self):                 
        self.linear1.bias.data.fill_(0)                                  # Initialize bias values at zero
        nn.init.xavier_uniform_(self.linear1.weight, gain=np.sqrt(2))    # Initialize weight values with xavier initialization


    # ----------------------------------------------------------------------------------------------------------------
    # Make a complete prediction (strains & derivatives)
    def predict(self): 

        if Function_flag == 1:
            data = self.X_TrainingData_SNS_Coarse_GP
        elif Function_flag == 0:
            data = self.X_PredictData_SNS_Coarse_GP    
        else:
            print("Check your Function_flag variable!") 

        predictions                         = self.forward(data)                                    # enonlocal 
        # enonlocal_grad_xy, enonlocal_elocal = self.first_partial_enonlocal(predictions, data)

        return predictions #, enonlocal_grad_xy, enonlocal_elocal


###########################################################################################################################################
# Initialize the model
model = Seq2Seq(Dimensions_input, Dimensions_output, PredictData_SNS_Coarse_GP_np, All_TrainingData_SNS_Coarse_GP_np)

if Function_flag == 1:

    # ----------------------------------------------------------------------------------------------------------------
    # Train the model
    Adam_start_time = time.time()
    model.train(Adam_epochs)
    total_training_time = time.time() - Adam_start_time
    Training_times_vec.append(total_training_time)

    # ----------------------------------------------------------------------------------------------------------------
    # Make a prediction 
    predictions = model.predict() #, enonlocal_grad_xy, enonlocal_elocal 
    
    # Un-scale predictions
    if Scaling_flag == "decimal":

        predictions         = torch.divide(predictions, decimal_scaling)
        # enonlocal_grad_xy   = torch.divide(enonlocal_grad_xy, decimal_scaling)
    
    elif Scaling_flag == "amsf":
        
        amsf_scaling = torch.unsqueeze(torch.tensor(amsf_scaling.astype(np.float32)), dim=-1).to(device)
        predictions         = torch.divide(predictions, amsf_scaling)
        # enonlocal_grad_xy   = torch.divide(enonlocal_grad_xy, amsf_scaling)

    elif Scaling_flag == "minmax":
    
        a = torch.unsqueeze(torch.tensor(a.astype(np.float32)), dim=-1).to(device)
        b = torch.unsqueeze(torch.tensor(b.astype(np.float32)), dim=-1).to(device)
        predictions       = torch.divide(predictions - b, a)
        # enonlocal_grad_xy = torch.divide(enonlocal_grad_xy, a)

    else:
        print("Check your Scaling_flag variable!")


    # ----------------------------------------------------------------------------------------------------------------
    # Print statements and saving
    print("Predictions shape: ", predictions.shape)
    print("Adam training time: ",   Training_times_vec[0], " sec")
    print("L-BFGS training time: ", Training_times_vec[1], " sec")
    print("Total training time: ",  Training_times_vec[2], " sec")
    print("N. L-BFGS iterations: ", Iter_LBFGS)

    torch.save(model.state_dict(), fileDir)

    scipy.io.savemat('TrainModeResults_' + Model_name + '.mat', 
                    {'predictions'              : predictions.detach().cpu().numpy(),
                    # 'enonlocal_grad_xy'         : enonlocal_grad_xy.detach().cpu().numpy(),
                    # 'enonlocal_elocal'          : enonlocal_elocal.detach().cpu().numpy(),                   
                    'training_times'            : Training_times_vec, 
                    'Enonlocal_term_norm_vec'   : Enonlocal_term_norm_vec,
                    'Laplacian_term_norm_vec'   : Laplacian_term_norm_vec,
                    'Elocal_term_norm_vec'      : Elocal_term_norm_vec,                                        
                    'loss_vec'                  : loss_vec, 
                    'loss_PDE_GP_vec'           : loss_PDE_GP_vec, 
                    'loss_PDEBCs_nodeslrb_vec'  : loss_PDEBCs_nodeslrb_vec, 
                    'loss_PDEBCs_nodesbtb_vec'  : loss_PDEBCs_nodesbtb_vec, 
                    'loss_Data_GP_vec'          : loss_Data_GP_vec,
                    'Iter_LBFGS'                : Iter_LBFGS})
    
elif Function_flag == 0:

    # ----------------------------------------------------------------------------------------------------------------
    # Load model and make a prediction on the increment of interest
    model.load_state_dict(torch.load(fileDir, map_location = device))

    predict_start_time = time.time()
    predictions = model.predict() # , enonlocal_grad_xy, enonlocal_elocal
    predict_end_time = time.time() - predict_start_time
    print(predict_end_time)
    
    # Un-scale predictions
    if Scaling_flag == "decimal":

        predictions         = torch.divide(predictions, decimal_scaling)
        # enonlocal_grad_xy   = torch.divide(enonlocal_grad_xy, decimal_scaling)

    elif Scaling_flag == "amsf":
        
        amsf_scaling = torch.unsqueeze(torch.tensor(amsf_scaling.astype(np.float32)), dim=-1).to(device)
        predictions         = torch.divide(predictions, amsf_scaling)
        # enonlocal_grad_xy   = torch.divide(enonlocal_grad_xy, amsf_scaling)

    elif Scaling_flag == "minmax":
    
        a = torch.unsqueeze(torch.tensor(a.astype(np.float32)), dim=-1).to(device)
        b = torch.unsqueeze(torch.tensor(b.astype(np.float32)), dim=-1).to(device)
        predictions       = torch.divide(predictions - b, a)
        # enonlocal_grad_xy = torch.divide(enonlocal_grad_xy, a)

    else:
        print("Check your Scaling_flag variable!")

    gc.collect()
    scipy.io.savemat('IFENN_Predictions_' + Model_name + '.mat', 
                    {'predictions'       : predictions.detach().cpu().numpy()})#,
                    #  'enonlocal_grad_xy' : enonlocal_grad_xy.detach().cpu().numpy(),
                    #  'enonlocal_elocal'  : enonlocal_elocal.detach().cpu().numpy()})

else:
    print("Check your Function_flag variable!") 

# def count_parameters(model):
#     return sum(p.numel() for p in model.parameters() if p.requires_grad)
# print(count_parameters(model))


        
