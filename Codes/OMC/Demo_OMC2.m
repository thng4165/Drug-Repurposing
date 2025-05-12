clear all
addpath('C:\Trang\KIProjects\ComprehensionDR\OMC_method\OMC-master');
tic;
rng('default')
myseed = 2024;
rng(myseed);

%% 1. Load Datesets

% load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\HDVD
% load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\LAGCN
% load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\Fdataset
% load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\Cdataset
% load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\LRSSL
% load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\Ydataset
%load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\DNdataset
% Wrr = drug;
% Wdd = disease;
% Wdr = didr;
% Wrd = Wdr';


%% CSV files
%%% oMat-MechDB dataset 

% disease = readmatrix('C:\Trang\KIProjects\ComprehensionDR\Datasets\CSVDatafiles\rare_disease_sim.csv');
% drug = readmatrix('C:\Trang\KIProjects\ComprehensionDR\Datasets\CSVDatafiles\rare_drug_sim.csv');
% didr = readmatrix('C:\Trang\KIProjects\ComprehensionDR\Datasets\CSVDatafiles\interact.csv');


%%% hsdn dataset
% disease = readmatrix('C:\Trang\KIProjects\ComprehensionDR\Datasets\CSVDatafiles\hsdn_disease_sim.csv');
% drug = readmatrix('C:\Trang\KIProjects\ComprehensionDR\Datasets\CSVDatafiles\hsdn_drug_sim.csv');
% didr = readmatrix('C:\Trang\KIProjects\ComprehensionDR\Datasets\CSVDatafiles\hsdn_interact.csv');
% 
% 
% Wrr = drug; % 150 drug
% Wdd = disease; % 89 disease
% Wdr = didr'; % 89*150
% Wrd = Wdr'; % 150*89


% a = load('C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\rare_disease_drug.mat');
% Wrd = a.data;
% Wdr = Wrd';
% b = load('C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\rare_drug_sim.mat');
% Wrr = b.data;
% c = load('C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\rare_disease_sim.mat');
% Wdd = c.data;



a = load('C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\hsdn_MechDB_dd_association_numeric.mat');
Wrd = a.data;
Wdr = Wrd';
b = load('C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\hsdn_MechDB_drug_sim.mat');
Wrr = b.data;
c = load('C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\hsdn_MechDB_disease_sim_GIP.mat');
Wdd = c.data;


%% Cross Validation


    nfolds = 10;
    positive_id1 = find(Wdr);
    positive_id2 = find(Wdr==0);
    crossval_idx1 = crossvalind('Kfold',positive_id1(:),nfolds);
    crossval_idx2 = crossvalind('Kfold',positive_id2(:),nfolds);
    crossval_idx = [crossval_idx1; crossval_idx2];

    % crossval_idx = crossvalind('Kfold', Wdr(:), nfolds); %% CV matrix disease*drug

%% 2. OMC2 algorithm
alpha = 1;
beta = 10;
K = 10;
tol1 = 2*1e-3;
tol2 = 1*1e-5;
maxiter = 300;
[nr,nc] = size(Wdr);


predicted_score = zeros(nr, nc); % final predicted matrix, fold by both 0 and 1
Null_dist = {};

for fold = 1:nfolds
    fold
    
    Tdat0.fold_number = fold;
    train_idx = find(crossval_idx ~= fold); 
    test_idx  = find(crossval_idx == fold);

    train_data = Wdr; 
    train_data(test_idx) = NaN;
    % Tdat0.train_set = train_data;
    % Tdat0.test_idx = test_idx;
    
    P_TMat = train_data;
    P_TMat(test_idx) = 0;


    row_no = find(sum(P_TMat, 2) == 0);

    
    
    if isempty(row_no) == 0
        P_TMat_new1 = KNN_diseaseS(P_TMat, Wdd, K);          %KNN Preprocessing
        P_TMat_new = P_TMat_new1 + P_TMat;
    else
        P_TMat_new =P_TMat;
    end

    T1 = [Wrr; P_TMat_new];
    [t1, t2] = size(T1);
    trIndex1 = double(T1 ~= 0);
    [W1, iter1] = BNNR(alpha, beta, T1, trIndex1, tol1, tol2, maxiter, 0, 1);
    M_ResultMat1 = W1((t1-nr+1):t1, 1:nc);

    col_no = find(sum(P_TMat, 1) == 0);
    if isempty(col_no) == 0
        P_TMat_new2 = KNN_drugS(P_TMat, Wrr, K);             %KNN Preprocessing
        P_TMat_new = P_TMat_new2 + P_TMat;
    else
        P_TMat_new = P_TMat;
    end
    T2 = [P_TMat_new, Wdd];
    [t_1, t_2] = size(T2);
    trIndex2 = double(T2 ~= 0);
    [W2, iter2] = BNNR(alpha, beta, T2, trIndex2, tol1, tol2, maxiter, 0, 1);
    M_ResultMat2 = W2(1:nr, 1:nc);

    M_recovery = (M_ResultMat1 + M_ResultMat2) / 2;     

    predicted_score(test_idx) = M_recovery(test_idx); 

    % Tdat0.train_approx = M_recovery;
    % Tdat0.test_approx = predicted_score(test_idx);
end
% writematrix(predicted_score, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\OMC_method\OMC-master\Figure1tun\OMC_HDVDdata_rerun.csv');
 % writematrix(predicted_score, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\OMC_method\OMC-master\Figure1tun\OMC_LAGCNdata_rerun.csv');
 % writematrix(predicted_score, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\OMC_method\OMC-master\Figure1tun\OMC_Fdata_rerun.csv');
 % writematrix(predicted_score, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\OMC_method\OMC-master\Figure1tun\OMC_Cdata_rerun.csv');
 % writematrix(predicted_score, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\OMC_method\OMC-master\Figure1tun\OMC_LRSSLdata_rerun.csv');
 % writematrix(predicted_score, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\OMC_method\OMC-master\Figure1tun\OMC_Ydata_rerun.csv');
 % writematrix(predicted_score, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\OMC_method\OMC-master\Figure1tun\OMC_oMatdata_rerun_neworder.csv');
 writematrix(predicted_score, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\OMC_method\OMC-master\Figure1tun\OMC_hsdndata_rerun_neworder.csv');
 % writematrix(predicted_score, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\OMC_method\OMC-master\Figure1tun\OMC_SCMFDDLdata_rerun.csv');
 %writematrix(predicted_score, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\OMC_method\OMC-master\Figure1tun\OMC_DNdata_rerun.csv');
 % writematrix(predicted_score, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\OMC_method\OMC-master\Figure1tun\OMC_iDrugdata_rerun.csv');
  %writematrix(predicted_score, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\OMC_method\OMC-master\Figure1tun\OMC_TLHGBIdata_rerun.csv');
%%  CV: AUC-AUPR
%%% sorting by disease
  % inputObs_matrix = didr'; % transpose to sort by column - disease
  % inputObs_matrix = didr; %for Omat and hsdn
  % prediction_matrix = predicted_score';

  %%% sorting by disease
  
  inputObs_matrix = Wrd; 
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
