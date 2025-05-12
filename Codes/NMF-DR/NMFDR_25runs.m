clear all;
addpath('C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\NMF-DR-main\Code');

tic;

%% Load Dataset

    
% load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\HDVD
% load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\LAGCN
% load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\Fdataset
% load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\Fdataset
%load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\Cdataset
% load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\LRSSL
  % load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\Ydataset
 % load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\iDrug

% load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\DNdataset
% load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\TLHGBI
 % load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\SCMFDD_L
    
    % spy(DiDrAMat)
    % title('Sparsity')
    % ylabel('Disease');
    % xlabel('Drugs');
   

% DrugSimMat = drug; 
% DiseaseSimMat = disease; 
% Wdr = didr; 
% DiDrAMat = Wdr';

%%
%% oMat-MechDB dataset 
% disease = readmatrix('C:\Trang\KIProjects\ComprehensionDR\Datasets\CSVDatafiles\rare_disease_sim.csv');
% drug = readmatrix('C:\Trang\KIProjects\ComprehensionDR\Datasets\CSVDatafiles\rare_drug_sim.csv');
% didr = readmatrix('C:\Trang\KIProjects\ComprehensionDR\Datasets\CSVDatafiles\interact.csv');

% %%% hsdn dataset
% disease = readmatrix('C:\Trang\KIProjects\ComprehensionDR\Datasets\CSVDatafiles\hsdn_disease_sim.csv');
% drug = readmatrix('C:\Trang\KIProjects\ComprehensionDR\Datasets\CSVDatafiles\hsdn_drug_sim.csv');
% didr = readmatrix('C:\Trang\KIProjects\ComprehensionDR\Datasets\CSVDatafiles\hsdn_interact.csv');
% % 
% 
% DrugSimMat = drug; 
% DiseaseSimMat = disease; 
% Wdr = didr; 
% DiDrAMat = Wdr;


a = load('C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\hsdn_MechDB_dd_association_numeric.mat');
Wdr = a.data;
DiDrAMat = Wdr;
b = load('C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\hsdn_MechDB_drug_sim.mat');
DrugSimMat = b.data;
c = load('C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\hsdn_MechDB_disease_sim_GIP.mat');
DiseaseSimMat = c.data;



    % a = load('C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\rare_disease_drug.mat');
    % Wdr = a.data;
    % DiDrAMat = Wdr;    
    % b = load('C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\rare_drug_sim.mat');
    % DrugSimMat = b.data;
    % c = load('C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\rare_disease_sim.mat');
    % DiseaseSimMat = c.data;


%%

r =50;
[m,n] = size(DiDrAMat); 
maxiter = 1e6; timelimit = 10;
nfolds = 10;


NMFDR_AUC = [];          
NMFDR_AUPR = [];    
seed_values = randsample(100000:200000, 2, false);

parfor iseed = 1:length(seed_values)

    
    fprintf('iseed: %d\n', iseed);
    myseed = seed_values(iseed);
    
    rng(myseed, 'twister');
    predicted_score = zeros(m,n); % final predicted matrix, fold by both 0 and 1
    % crossval_idx = crossvalind('Kfold', DiDrAMat(:), nfolds); %% CV matrix disease*drug

    positive_id1 = find(DiDrAMat);
    positive_id2 = find(DiDrAMat==0);
    crossval_idx1 = crossvalind('Kfold',positive_id1(:),nfolds);
    crossval_idx2 = crossvalind('Kfold',positive_id2(:),nfolds);
    crossval_idx = [crossval_idx1; crossval_idx2];

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
         % [W2,H2] = SVDNMF(X,r); 
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
        % [W2n,H2n,e2n] = HALSacc(X,W2,H2,0,0,maxiter,timelimit); 
         [W3n,H3n,e3n] = HALSacc(X,W3,H3,0,0,maxiter,timelimit); 
        %[W4n,H4n,e4n] = HALSacc(X,W4,H4,0,0,maxiter,timelimit); 
        
        %%-------------------
         % Y2=W2n*H2n;
         Y3=W3n*H3n;
        %Y2=W4n*H4n;
        % predicted_score(test_idx) = Y2(test_idx); 
        predicted_score(test_idx) = Y3(test_idx); 
        end


    %%  CV: AUC-AUPR
    %%% sorting by disease
      inputObs_matrix = DiDrAMat ; % transpose to sort by column - disease
     prediction_matrix = predicted_score;
    
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
    
    
    rr = size(inputObs_matrix, 1);
    
    %% Calculate AUC_AUPR
    
            for cutoff=1:size(inputObs_matrix, 1)
                P_matrix = sorted_inputObs_matrix(1:cutoff, :);
                %N_matrix = sorted_inputObs_matrix((cutoff+1):rr, :);
                if cutoff < rr
                    N_matrix = sorted_inputObs_matrix((cutoff+1):rr, :);
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
            end
    
        % Compute AUC and AUPR
        AUC = trapz(fpr_list, tpr_list); % AUC
        AUPR = trapz(recall_list, precision_list); % AUPR
        NMFDR_AUC(iseed) = AUC;
        NMFDR_AUPR(iseed) = AUPR;
end

% writematrix(NMFDR_AUC, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\NMF-DR-main\Code\Result25runs\NMFDR_SCMFDD_L_AUC_25runs_rerun.csv');
% writematrix(NMFDR_AUPR, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\NMF-DR-main\Code\Result25runs\NMFDR_SCMFDD_L_AUPR_25runs_rerun.csv');

% writematrix(NMFDR_AUC, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\NMF-DR-main\Code\Result25runs\NMFDR_LAGCN_AUC_25runs_rerun.csv');
% writematrix(NMFDR_AUPR, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\NMF-DR-main\Code\Result25runs\NMFDR_LAGCN_AUPR_25runs_rerun.csv');

% writematrix(NMFDR_AUC, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\NMF-DR-main\Code\Result25runs\NMFDR_Fdata_AUC_25runs_rerun.csv');
% writematrix(NMFDR_AUPR, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\NMF-DR-main\Code\Result25runs\NMFDR_Fdata_AUPR_25runs_rerun.csv');

% writematrix(NMFDR_AUC, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\NMF-DR-main\Code\Result25runs\NMFDR_Cdata_AUC_25runs_rerun.csv');
% writematrix(NMFDR_AUPR, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\NMF-DR-main\Code\Result25runs\NMFDR_Cdata_AUPR_25runs_rerun.csv');
% writematrix(NMFDR_AUC, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\NMF-DR-main\Code\Result25runs\NMFDR_Ydata_AUC_25runs_rerun.csv');
% writematrix(NMFDR_AUPR, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\NMF-DR-main\Code\Result25runs\NMFDR_Ydata_AUPR_25runs_rerun.csv');


% writematrix(NMFDR_AUC, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\NMF-DR-main\Code\Result25runs\NMFDR_hsdn_AUC_25runs_rerun.csv');
% writematrix(NMFDR_AUPR, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\NMF-DR-main\Code\Result25runs\NMFDR_hsdn_AUPR_25runs_rerun.csv');


writematrix(NMFDR_AUC, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\NMF-DR-main\Code\Result25runs\NMFDR_hsdn_AUC_2runs_rerun_neworder.csv');
writematrix(NMFDR_AUPR, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\NMF-DR-main\Code\Result25runs\NMFDR_hsdn_AUPR_2runs_rerun_neworder.csv');

t = toc;