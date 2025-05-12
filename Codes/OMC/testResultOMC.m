clear all
addpath('C:\Trang\KIProjects\ComprehensionDR\OMC_method\OMC-master');

%% 1. Load Datesets
%load Fdataset
%load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\HDVD 

%load LAGCN
%load Fdataset
%load Cdataset
%load LRSSL
%load Ydataset
%load DNdataset
load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\SCMFDD_L
%load iDrug
%load TLHGBI

Wrr = drug;
Wdd = disease;
Wdr = didr;
Wrd = Wdr';


%% CSV files
%%% oMat-MechDB dataset 
% disease = readmatrix('rare_disease_sim.csv');
% drug = readmatrix('rare_drug_sim.csv');
% didr = readmatrix('interact.csv');

%%% hsdn dataset
% disease = readmatrix('hsdn_disease_sim.csv');
% drug = readmatrix('hsdn_drug_sim.csv');
% didr = readmatrix('hsdn_interact.csv');


% Wrr = drug; % 150 drug
% Wdd = disease; % 89 disease
% Wdr = didr'; % 89*150
% Wrd = Wdr'; % 150*89
% 




    nfolds = 10;
    myseed = 2024;
    rng(myseed, 'twister');

 % readmatrix('OMC_HDVDdata.csv');
 % readmatrix('OMC_LAGCNdata.csv');
 % readmatrix('OMC_Fdata.csv');
 % readmatrix('OMC_Cdata.csv');
 % readmatrix('OMC_LRSSLdata.csv');
 % readmatrix('OMC_Ydata.csv');
 % readmatrix('OMC_oMatdata.csv');
 % readmatrix('OMC_hsdndata.csv');
 predicted_score = readmatrix('C:\Trang\KIProjects\ComprehensionDR\OMC_method\OMC-master\1runs\OMC_SCMFDDLdata.csv');
%  readmatrix('OMC_DNdata.csv');
%  readmatrix('OMC_iDrugdata.csv');
%  readmatrix('OMC_TLHGBIdata.csv');
% %%  CV: AUC-AUPR
%%% sorting by disease
  % inputObs_matrix = didr'; % transpose to sort by column - disease
  % inputObs_matrix = didr; %for Omat and hsdn
  % prediction_matrix = predicted_score';

  %%% sorting by disease
  inputObs_matrix = didr'; % transpose to sort by column - disease
  % inputObs_matrix = didr; %for Omat and hsdn
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


