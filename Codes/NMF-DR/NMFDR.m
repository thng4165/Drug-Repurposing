clear all;
addpath('NMF-DR');

tic;


%% Load Dataset
 load Datasets\MatlabDataFiles\HDVD

DrugSimMat = drug; 
DiseaseSimMat = disease; 
Wdr = didr; 
DiDrAMat = Wdr';

% oMat-MechDB dataset 
% a = load('C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\hsdn_MechDB_dd_association_numeric.mat');
% Wdr = a.data;
% DiDrAMat = Wdr;
% b = load('C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\hsdn_MechDB_drug_sim.mat');
% DrugSimMat = b.data;
% c = load('C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\hsdn_MechDB_disease_sim_GIP.mat');
% DiseaseSimMat = c.data;

% hsdn dataset
 a = load('C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\rare_disease_drug.mat');
 Wdr = a.data;
 DiDrAMat = Wdr;    
 b = load('C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\rare_drug_sim.mat');
 DrugSimMat = b.data;
 c = load('C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\rare_disease_sim.mat');
 DiseaseSimMat = c.data;


%% NMF-DR algorithm
% Parameters
r =50;
[m,n] = size(DiDrAMat); 
maxiter = 1e6; timelimit = 10;
nfolds = 10;

myseed = 2024;
rng(myseed);

% Cross Validation
 positive_id1 = find(DiDrAMat);
 positive_id2 = find(DiDrAMat==0);
 crossval_idx1 = crossvalind('Kfold',positive_id1(:),nfolds);
 crossval_idx2 = crossvalind('Kfold',positive_id2(:),nfolds);
 crossval_idx = [crossval_idx1; crossval_idx2];
  
 % NMFDR
 predicted_score = zeros(m,n); % final predicted matrix, fold by both 0 and 1

for fold = 1:nfolds
    fold
    
    train_idx = find(crossval_idx ~= fold); 
    test_idx  = find(crossval_idx == fold);

    train_data = DiDrAMat; % train matrix
    train_data(test_idx) = 0;
   
    X=DrugSimMat* train_data* DiseaseSimMat;
    X = train_data .* X;
    X = X/max(X(:));

% Compare three SVD-based NMF initializations and random initialization 
% % 1) NNDSVD
% fprintf('Running NNDSVD...'), 
% tic; 
% [W1,H1] = NNDSVD(X,r,1);
% fprintf(' Done. Computational time = %2.2f s.\n', toc); 
% % 2) SVD-NMF
% fprintf('Running SVD-NMF...'), 
% tic; 
 %[W2,H2] = SVDNMF(X,r); 
% fprintf(' Done. Computational time = %2.2f s.\n', toc); 
% % % 3) NNSVDLRC
% fprintf('Running NNSVD-LRC...'), 
% tic; 
 [W3,H3] = NNSVDLRC(X,r); 
% fprintf(' Done. Computational time = %2.2f s.\n', toc); 
% % 4) Random initialization
% W4 = rand(m,r); 
% H4 = rand(r,n); 

% Display evolution of A-HALS error for the 3 initializations 
% [W1n,H1n,e1n] = HALSacc(X,W1,H1,0,0,maxiter,timelimit); 
%[W2n,H2n,e2n] = HALSacc(X,W2,H2,0,0,maxiter,timelimit); 
  [W3n,H3n,e3n] = HALSacc(X,W3,H3,0,0,maxiter,timelimit); 
%[W4n,H4n,e4n] = HALSacc(X,W4,H4,0,0,maxiter,timelimit); 

%%-------------------
  %Y2=W2n*H2n;
 Y3=W3n*H3n;
%Y2=W4n*H4n;
% predicted_score(test_idx) = Y2(test_idx); 

 predicted_score(test_idx) = Y3(test_idx); 
 
end

%%  CV: AUC-AUPR
% sorting by disease
 inputObs_matrix = DiDrAMat ; % transpose to sort by column - disease
 prediction_matrix = predicted_score;

% sort inputObs_matrix by column using the decreasing order by column of prediction_matrix
  res = sort_matrix(prediction_matrix, inputObs_matrix);
  sorted_inputObs_matrix = res.y_sorted;
  sorted_score_matrix = res.score_sorted;
  sort_index = res.sort_index;

 % Initialize lists
  tpr_list = [];
  fpr_list = [];
  recall_list = [];
  precision_list = [];


rr = size(inputObs_matrix, 1);

% Calculate AUC_AUPR

        for cutoff=1:size(inputObs_matrix, 1)
            P_matrix = sorted_inputObs_matrix(1:cutoff, :);
            %N_matrix = sorted_inputObs_matrix((cutoff+1):rr, :);
            if cutoff < rr
                N_matrix = sorted_inputObs_matrix((cutoff+1):rr, :);
            else
                N_matrix = [];  % If cutoff = rr, set N_matrix to empty
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
