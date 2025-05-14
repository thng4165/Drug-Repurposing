clear all;
rng('default'); seed = 2024; rng(seed)
addpath("HGIMC");

addpath('Functions');
tic;
%% 1. Load Datesets

load Datasets\MatlabDataFiles\HDVD
% load Datasets\MatlabDataFiles\LAGCN
% load Datasets\MatlabDataFiles\Fdataset
% load Datasets\MatlabDataFiles\Cdataset.mat
% load Datasets\MatlabDataFiles\Cdataset
% load Datasets\MatlabDataFiles\LRSSL
% load Datasets\MatlabDataFiles\Ydataset

 R = drug; 
 D = disease; 
 A_DR_original = didr; 


% oMat-MechDB dataset
% a = load('Datasets\MatlabDataFiles\rare_disease_drug.mat');
% didr = a.data;
% A_DR_original = didr';
% b = load('Datasets\MatlabDataFiles\rare_drug_sim.mat');
% R = b.data;
% c = load('Datasets\MatlabDataFiles\rare_disease_sim.mat');
% D = c.data;

% hsdn-MechDB dataset
% a = load('Datasets\MatlabDataFiles\hsdn_MechDB_dd_association_numeric.mat');
% didr = a.data;
% A_DR_original = didr';
% b = load('Datasets\MatlabDataFiles\hsdn_MechDB_drug_sim.mat');
% R = b.data;
% c = load('Datasets\MatlabDataFiles\hsdn_MechDB_disease_sim_GIP.mat');
% D = c.data;


%% Checking with HGIMC datasets
% load Datasets\HGIMCdata\Fdataset_ms
% load Datasets\HGIMCdata\Cdataset_ms
% load Datasets\HGIMCdata\Ydataset_ms
% load Datasets\HGIMCdata\iDrug_ms

% A_DR_original = didr;
% R = (drug_ChemS+drug_AtcS+drug_SideS+drug_DDIS+drug_TargetS)/5;
% D = (disease_PhS+disease_DoS)/2;



%% Cross Validation

%  10 folds CV for both 0 (unknown asso) and 1 (known asso)

    nfolds = 10;
    positive_id1 = find(A_DR_original);
     positive_id2 = find(A_DR_original==0);
     crossval_idx1 = crossvalind('Kfold',positive_id1(:),nfolds);
     crossval_idx2 = crossvalind('Kfold',positive_id2(:),nfolds);
     crossval_idx = [crossval_idx1; crossval_idx2];
    


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

 predicted_score = zeros(nr, nc); % final predicted matrix, fold by both 0 and 1
 


 for fold = 1:nfolds
    fold
    

    train_idx = find(crossval_idx ~= fold); 
    test_idx  = find(crossval_idx == fold);

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

    predicted_score(test_idx) = A_recovery(test_idx); 

    
 end

%%  CV: AUC-AUPR
%%% sorting by disease
  inputObs_matrix = A_DR_original';
  prediction_matrix = predicted_score';

%%% sort inputObs_matrix by column using the decreasing order by column of prediction_matrix
  res = sort_matrix(prediction_matrix, inputObs_matrix);
  sorted_inputObs_matrix = res.y_sorted;
  sorted_score_matrix = res.score_sorted;
  sort_index = res.sort_index;

 %%% Initialize lists
  tpr_list = [];
  fpr_list = [];
  recall_list = [];
  precision_list = [];


r = size(inputObs_matrix, 1);

%% Calculate AUC_AUPR

        for cutoff=1:size(inputObs_matrix, 1)
            P_matrix = sorted_inputObs_matrix(1:cutoff, :);
            %N_matrix = sorted_inputObs_matrix((cutoff+1):r, :);
            if cutoff < r
                N_matrix = sorted_inputObs_matrix((cutoff+1):r, :);
            else
                N_matrix = [];  % If cutoff = r, set N_matrix to empty
            end



            TP = sum(P_matrix(:) == 1);
            FP = sum(P_matrix(:) == 0);
            TN = sum(N_matrix(:) == 0);
            FN = sum(N_matrix(:) == 1);
            tpr = TP / (TP + FN);
            fpr = FP / (FP + TN);
            recall_ = TP / (TP + FN);
            precision_ = TP / (TP + FP);
            accuracy_ = (TN + TP) / (TN + TP + FN + FP);
            f1_ = (2 * TP) / (2 * TP + FP + FN);
            tpr_list = [tpr_list,tpr];
            fpr_list = [fpr_list,fpr];
            recall_list = [recall_list,recall_];
            precision_list = [precision_list,precision_];
           % accuracy_list = [accuracy_list,accuracy_];
           % F1_list = [F1_list,f1_];
        end



    % Compute AUC and AUPR
    AUC = trapz(fpr_list, tpr_list); % AUC
    AUPR = trapz(recall_list, precision_list); % AUPR

    
    figure;
    subplot(1, 2, 1);
    plot(fpr_list, tpr_list, 'LineWidth', 2);
    title(['AUC = ', num2str(AUC)], 'FontSize', 14);
    xlabel('FPR (1-specificity)', 'FontSize', 12);
    ylabel('TPR (sensitivity)', 'FontSize', 12);
    grid on;
    hold on;
    plot([0, 1], [0, 1], '--', 'Color', [0.5, 0.5, 0.5]);
    hold off;

    subplot(1, 2, 2);
    plot(recall_list, precision_list, 'LineWidth', 2);
    title(['AUPR = ', num2str(AUPR)], 'FontSize', 14);
    xlabel('Recall', 'FontSize', 12);
    ylabel('Precision', 'FontSize', 12);
    grid on;

t = toc;
