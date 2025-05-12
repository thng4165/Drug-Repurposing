clear all;
rng('default'); seed = 12345; rng(seed)
addpath("C:\Trang\KIProjects\ComprehensionDR\HGIMC");

addpath('Functions');
%% 1. Load Datesets
 % load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\Cdataset.mat
%load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\HDVD
% load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\LAGCN
 load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\Fdataset
    % load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\Cdataset
     % load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\LRSSL
    % load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\Ydataset
   % load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\DNdataset
    % load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\SCMFDD_L

    %load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\iDrug



R = drug; 
D = disease; 
A_DR_original = didr; 




%% Cross Validation

% % 10 folds CV for both 0 (unknown asso) and 1 (known asso)
% % almost same elements in each folds for both 0 and 1
    nfolds = 10;
 
%% 2. HGIMC algorithm
%% Parameter
[nr,nc] = size(A_DR_original);

alpha = 10; 
beta = 10; 
gamma = 0.1; 
threshold = 0.1;
maxiter = 300; 
tol1 = 2*1e-3;   
tol2 = 1*1e-5;

%% HGIMC
[positiveId, crossval_id] = train_test_split(A_DR_original, nfolds, '1');
    O = [];
    P = [];
  for fold = 1:nfolds
      fold
        PtestID = positiveId(crossval_id == fold);
        negativeID = find(A_DR_original == 0);
        num = numel(negativeID);
        Nidx = randperm(num);
        NtestID = negativeID(Nidx(1:length(PtestID)));
           
    test_idx  = PtestID;

    train_data = A_DR_original; 
    train_data(test_idx) = NaN;
       
    A_DR = train_data; 
    A_DR(test_idx) = 0; 
    
    % 2.1 Bounded Matrix Completion
    trIndex = double(A_DR ~= 0);
    [A_bmc, iter] = fBMC(alpha, beta, A_DR, trIndex, tol1, tol2, maxiter, 0, 1);
    A_DR0 = A_bmc.*double(A_bmc > threshold);

    % 2.2 Gaussian Radial Basis function
    A_RR = fGRB(R, 0.5);
    A_DD = fGRB(D, 0.5);

    % 2.3 Heterogeneous Graph Inference 
    A_recovery = fHGI(gamma, A_DD, A_RR, A_DR0);

    origin = [A_DR_original(PtestID); A_DR_original(NtestID)];
    pred = [A_recovery(PtestID); A_recovery(NtestID)];
    O = [O;origin];
    P = [P;pred];
      
  end

 inputObs_matrix = O; % transpose to sort by column - disease
 prediction_matrix = P;



[roc_x, roc_y, ~, AUC] = perfcurve(inputObs_matrix, prediction_matrix, true);

[pr_x, pr_y, ~, AUPR] = perfcurve(inputObs_matrix, prediction_matrix, 1, 'xCrit', 'reca', 'yCrit', 'prec');

writematrix(roc_x, 'HGIMC_roc_x_Fdata.csv');
writematrix(roc_y, 'HGIMC_roc_y_Fdata.csv');
writematrix(pr_x, 'HGIMC_pr_x_Fdata.csv');
writematrix(pr_y, 'HGIMC_pr_y_Fdata.csv');

% Plot ROC curve
figure;
plot(roc_x, roc_y, 'b-', 'LineWidth', 2);
xlabel('False Positive Rate','FontSize', 12);
ylabel('True Positive Rate','FontSize', 12);
title(['ROC Curve (AUC = ' num2str(AUC) ')'], 'FontSize', 14);
grid on;
hold on;
plot([0, 1], [0, 1], '--', 'Color', [0.5, 0.5, 0.5]);
hold off;


figure;
plot(pr_x, pr_y, 'b-', 'LineWidth', 2);
xlabel('Recall', 'FontSize', 12);
ylabel('Precision', 'FontSize', 12);
title(['Precision-Recall Curve (AUPR = ' num2str(AUPR) ')'], 'FontSize', 14);
grid on;

