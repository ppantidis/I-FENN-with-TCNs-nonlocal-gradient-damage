This repository contains the code and data used in this article: 

Panos Pantidis, Habiba Eldababy, Diab Abueidda, Mostafa E. Mobasher (2024), "I-FENN with Temporal Convolutional Networks: expediting the load-history analysis of non-local gradient damage propagation"
https://arxiv.org/abs/2402.05460

The main MATLAB script is "FEM_Main_Script.m", from where either an I-FENN or FEM analysis can be executed.

The ".txt" files contain the information for the examples presented in the paper. The mesh details (nodes, elements, BCs) are located in the "Modelname_mesh.txt" files, and the analysis parameters are found in the "Modelname_parameters.txt").

The "TCN_IFENN.py" file has the code for the Temporal Convolutional Network, and it is used for both the offline training and its predictions within I-FENN. 



