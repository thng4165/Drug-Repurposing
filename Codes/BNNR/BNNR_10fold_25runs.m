clear all;
addpath('BNNR');


%parpool(24)

t0 = tic;
%% 1. Load Datesets
load Datasets\MatlabDataFiles\HDVD

Wrr = drug; 
Wdd = disease; 
Wdr = didr; 
Wrd = Wdr';

% oMat - MechDB dataset
% a = load('C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\rare_disease_drug.mat');
% Wrd = a.data;
% Wdr = Wrd';
% b = load('C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\rare_drug_sim.mat');
% Wrr = b.data;
% c = load('C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\rare_disease_sim.mat');
% Wdd = c.data;

% hsdn-MechDB dataset
% a = load('C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\hsdn_MechDB_dd_association_numeric.mat');
% Wrd = a.data;
% Wdr = Wrd';
% b = load('C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\hsdn_MechDB_drug_sim.mat');
% Wrr = b.data;
% c = load('C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\hsdn_MechDB_disease_sim_GIP.mat');
% Wdd = c.data;


%% Parameter

[dn,dr] = size(Wdr);

% parameter to get the good results.
maxiter = 300;
alpha = 1; 
beta = 10; 
tol1 = 2*1e-3;
tol2 = 1*1e-5;
nfolds = 10;

% Initialize
BNNR_AUC = [];           
BNNR_AUPR = [];          

% Generate 25 unique random integers
seed_values = randsample(100000:200000, 25, false);  


parfor iseed = 1:length(seed_values)
    fprintf('iseed: %d\n', iseed);
    
%% Cross Validation

% 10 folds CV for both 0 (unknown asso) and 1 (known asso)
 
    myseed = seed_values(iseed);
    rng(myseed,'twister');

    positive_id1 = find(Wdr);
    positive_id2 = find(Wdr==0);
    crossval_idx1 = crossvalind('Kfold',positive_id1(:),nfolds);
    crossval_idx2 = crossvalind('Kfold',positive_id2(:),nfolds);
    crossval_idx = [crossval_idx1; crossval_idx2];


%% BNNR models
    predicted_score = zeros(size(Wdr,1), size(Wdr,2)); % final predicted matrix, fold by both 0 and 1
    theta = 0;
for fold = 1:nfolds
    fold
    
    train_idx = find(crossval_idx ~= fold); 
    test_idx  = find(crossval_idx == fold);

    train_data = Wdr; % train matrix
    train_data(test_idx) = NaN;
    
  
    WRdr = train_data; 
    WRdr(test_idx) = theta;
    
    T = [Wrr, WRdr'; WRdr, Wdd];
    [t1, t2] = size(T);
    trIndex = double(T ~= 0);
    [WW,iter] = BNNR(alpha, beta, T, trIndex, tol1, tol2, maxiter, 0, 1);
    M_recovery = WW((t1-dn+1) : t1, 1 : dr);
    %M_recovery = (M_recovery');
    predicted_score(test_idx) = M_recovery(test_idx); 
    
end


%%  CV: AUC-AUPR
% sorting by disease
 inputObs_matrix = Wdr'; % transpose to sort by column - disease
 prediction_matrix = predicted_score';

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
        end

    % Compute AUC and AUPR
    AUC = trapz(fpr_list, tpr_list); % AUC
    AUPR = trapz(recall_list, precision_list); % AUPR
    BNNR_AUC(iseed) = AUC;
    BNNR_AUPR(iseed) = AUPR;
        
 end 

t1 = toc(t0);