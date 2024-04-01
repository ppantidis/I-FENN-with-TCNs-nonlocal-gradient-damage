clc
clear all
close all
format compact

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load TrainingData_StrB100_GP_time_features_lf
load sifenn

Pred1 = TrainingData_StrB100_GP_time_features_lf(:,1:41,:);
Pred2 = TrainingData_StrB100_GP_time_features_lf(:,26,:);

Pred3 = TrainingData_StrB100_GP_time_features_lf(:,1:26,:);
Pred3(:,1:25,:) = NaN;

Pred4 = TrainingData_StrB100_GP_time_features_lf(:,1:41,:);
Pred4(:,1:25,:)  = NaN;
Pred4(:,27:41,:) = NaN;

Pred5 = TrainingData_StrB100_GP_time_features_lf(:,1:41,:);
Pred5(:,27:41,:) = NaN;

Pred6 = TrainingData_StrB100_GP_time_features_lf(:,1:41,:);
Pred6(:,1,:) = NaN;

sifennPred1 = sifenn(1,1:41);
sifennPred2 = sifenn(1,26);
sifennPred3 = sifenn(1,1:26);
sifennPred4 = sifenn(1,1:41);
sifennPred5 = sifenn(1,1:41);
sifennPred6 = sifenn(1,1:41);


% -------------------------------------------------------------------------
Data_CP1 = TrainingData_StrB100_GP_time_features_lf(:,1:41,:);
Data_CP1(:,12:41,:) = NaN;
sifennCP1 = sifenn(1,1:41);

Data_CP2 = TrainingData_StrB100_GP_time_features_lf(:,1:41,:);
Data_CP2(:,22:41,:) = NaN;
sifennCP2 = sifenn(1,1:41);

Data_CP3 = TrainingData_StrB100_GP_time_features_lf(:,1:41,:);
Data_CP3(:,1:5,:) = NaN;
Data_CP3(:,12:41,:) = NaN;
sifennCP3 = sifenn(1,1:41);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save Pred1.mat Pred1
save Pred2.mat Pred2
save Pred3.mat Pred3
save Pred4.mat Pred4
save Pred5.mat Pred5
save Pred6.mat Pred6
save Data_CP1.mat Data_CP1
save Data_CP2.mat Data_CP2
save Data_CP3.mat Data_CP3

save sifennPred1.mat sifennPred1
save sifennPred2.mat sifennPred2
save sifennPred3.mat sifennPred3
save sifennPred4.mat sifennPred4
save sifennPred5.mat sifennPred5
save sifennPred6.mat sifennPred6
save sifennCP1.mat sifennCP1
save sifennCP2.mat sifennCP2
save sifennCP3.mat sifennCP3
