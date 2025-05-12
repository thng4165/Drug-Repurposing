clear all

% rng('default'); seed = 12345; rng(seed)
tic;

methodset  = {'VDA-GMSBMF'};

    %nCV = 5; 
    nfold  = 10; CVtype = 'CVa';    
  
    %% 1. Load Datesets
     % load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\HDVD    
    % load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\LAGCN
     % load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\Fdataset
    % load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\Cdataset
   %load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\LRSSL
    % load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\Ydataset
    % load DNdataset
    % load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\SCMFDD_L
    % load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\iDrug
    %load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\DNdataset

    % Wdv = didr'; %load([datadir,filesep,'matDrugVirus.txt']   ); 
    % Wvv = disease; %load([datadir,filesep,'matVirusVirus.txt']   ); 
    % Wdd = drug; %%load([datadir,filesep,'matDrugDrug.txt']  ); 

    %% oMat-MechDB dataset 
    % disease = readmatrix('C:\Trang\KIProjects\ComprehensionDR\Datasets\CSVDatafiles\rare_disease_sim.csv');
    % drug = readmatrix('C:\Trang\KIProjects\ComprehensionDR\Datasets\CSVDatafiles\rare_drug_sim.csv');
    % didr = readmatrix('C:\Trang\KIProjects\ComprehensionDR\Datasets\CSVDatafiles\interact.csv');

    % drug = readmatrix('C:\Trang\KIProjects\ComprehensionDR\Datasets\CSVDatafiles\hsdn_drug_sim.csv');
    % disease = readmatrix('C:\Trang\KIProjects\ComprehensionDR\Datasets\CSVDatafiles\hsdn_disease_sim.csv');
    % didr = readmatrix('C:\Trang\KIProjects\ComprehensionDR\Datasets\CSVDatafiles\hsdn_interact.csv');
    % 
    % 
    % Wdv = didr; 
    % Wvv = disease; 
    % Wdd = drug; 

 


    a = load('C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\hsdn_MechDB_dd_association_numeric.mat');
    Wdv = a.data;
    
    b = load('C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\hsdn_MechDB_drug_sim.mat');
    Wdd = b.data;
    c = load('C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\hsdn_MechDB_disease_sim_GIP.mat');
    Wvv = c.data;

    [Nd,Nv] = size(Wdv); 
    
VDA_AUC = [];
VDA_AUPR = [];
seed_values = randsample(100000:200000, 25, false); 
parfor iseed = 1:length(seed_values)
 
    fprintf('iseed: %d\n', iseed);
MatPredict = zeros(Nd,Nv);
%% Cross Validation

    nfolds = 10;
    seed = seed_values(iseed);
    % seed = 2024;
    rng('default');
rng(seed);
    for i_mtd = 1:length( methodset )
        method      = methodset{i_mtd} ;
        MatVecSet   = []; 
        tbScalarSet = [];   
        ii_fold = 0; 
		
        i_cv = 1;
            %disp([datadir,'--',method,':cv:',num2str(i_cv)]) 
            MatPredict=zeros( size( Wdv ) );
            % 
            [CVdata] = getKfoldCrossValidMatIndSet(Wdv, nfold,CVtype,'Unlabel', seed+i_cv ) ;         
            %     [CVdata] = getKfoldCrossValidMatIndSet(Wdv, nfold,CVtype,'Unlabel'  ) ;           
            IndSet_pos_test  = CVdata.MatIndSet_pos_test ; 
            IndSet_neg_test  = CVdata.MatIndSet_neg_test ; 
            IndSet_pos_train = CVdata.MatIndSet_pos_train; 
            IndSet_neg_train = CVdata.MatIndSet_neg_train ;                  
            % % % % % % % % % % % % 
            for i_fold=1: nfold
                i_fold
				%disp([datadir,'--',method,': i_cv-',num2str(i_cv),'; i_fold-',num2str(i_fold)])
                Ind_pos_test = IndSet_pos_test{i_fold} ;   
                Ind_neg_test = IndSet_neg_test{i_fold} ; 
                Ind_test     = union(Ind_pos_test,Ind_neg_test) ; %%% [Ind_neg_test;Ind_pos_test];
                matDV           = Wdv;
                matDV(Ind_test) = 0; 
                %
                Sdd=Wdd;
                Svv=Wvv;
                M_recovery = []; 
                % 
                switch method
                    case 'VDA-MSBMF'    %%% OK
                        lambda1 = 0.1; lambda2 = 0.1; lambda3 = lambda2; k = floor(min(Nd, Nv) * 0.7); tol1 = 2*1e-3;tol2 = 1*1e-4;      
                        [U, V, iter] = A_MSBMF(matDV, Sdd, Svv, lambda1, lambda2, lambda3, k, tol1, tol2, 300);
                        M_recovery = U * V';    

                    case 'VDA-GMSBMF'  %%% OK
                        %%%lambda1 = 0.1;lambda2 = 0.1;lambda3 = lambda2;k = floor(min_mn * 0.7); maxiter = 300;tol1 = 2*1e-3;tol2 = 1*1e-4;w=0.5;
                        %switch datadir
                         %   case 'VDdataset1'; gm= 0.5; w=0.3;lambda1 = 1;     lambda2 = lambda1; lambda3 = lambda2;  tol1 = 2*1e-30;tol2 = 1*1e-40;   
                         %   otherwise 
                                % gm= 0.5; w=0.4;lambda1 = 0.5; lambda2 =
                                % lambda1; lambda3 = lambda2;  tol1 =
                                % 2*1e-30;tol2 = 1*1e-40; % first run
                                gm= 0.5; w=0.4;lambda1 = 1; lambda2 = lambda1; lambda3 = lambda2;  tol1 = 2*1e-30;tol2 = 1*1e-40; %rerun
                        %end                          
                        [M_recovery] = A_VDA_GMSBMF(matDV, Sdd, Svv, gm,w, lambda1, lambda2, lambda3, floor( min(size(matDV)) * 1), tol1, tol2, 400) ; 
                            
                    otherwise; error(['There is no defintion: ', method]);                  
                end 

                % 
                MatPredict(Ind_test)= M_recovery(Ind_test);   
                
            end 

       

    end 
    predicted_score = MatPredict;

%%  CV: AUC-AUPR
%%% sorting by disease
  inputObs_matrix = Wdv; % transpose to sort by column - disease
 prediction_matrix = MatPredict;



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


    % figure;
    % subplot(1, 2, 1);
    % plot(fpr_list, tpr_list, 'LineWidth', 2);
    % title(['AUC = ', num2str(AUC)], 'FontSize', 14);
    % xlabel('FPR (1-specificity)', 'FontSize', 12);
    % ylabel('TPR (sensitivity)', 'FontSize', 12);
    % grid on;
    % hold on;
    % plot([0, 1], [0, 1], '--', 'Color', [0.5, 0.5, 0.5]);
    % hold off;
    % 
    % subplot(1, 2, 2);
    % plot(recall_list, precision_list, 'LineWidth', 2);
    % title(['AUPR = ', num2str(AUPR)], 'FontSize', 14);
    % xlabel('Recall', 'FontSize', 12);
    % ylabel('Precision', 'FontSize', 12);
    % grid on;




VDA_AUC(iseed) = AUC;
VDA_AUPR(iseed) = AUPR;

end
   

% writematrix(VDA_AUC, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\VDA_GMSBMF\VDA_GMSBMF-main\VDA_25runs\VDA_HDVD_AUC_25runs_rerun.csv');
% writematrix(VDA_AUPR, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\VDA_GMSBMF\VDA_GMSBMF-main\VDA_25runs\VDA_HDVD_AUPR_25runs_rerun.csv');

%writematrix(VDA_AUC, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\VDA_GMSBMF\VDA_GMSBMF-main\VDA_25runs\VDA_LAGCN_AUC_25runs_rerun.csv');
%writematrix(VDA_AUPR, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\VDA_GMSBMF\VDA_GMSBMF-main\VDA_25runs\VDA_LAGCN_AUPR_25runs_rerun.csv');

% writematrix(VDA_AUC, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\VDA_GMSBMF\VDA_GMSBMF-main\VDA_25runs\VDA_Fdata_AUC_25runs_rerun.csv');
% writematrix(VDA_AUPR, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\VDA_GMSBMF\VDA_GMSBMF-main\VDA_25runs\VDA_Fdata_AUPR_25runs_rerun.csv');
% writematrix(VDA_AUC, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\VDA_GMSBMF\VDA_GMSBMF-main\VDA_25runs\VDA_oMat_AUC_25runs_rerun.csv');
% writematrix(VDA_AUPR, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\VDA_GMSBMF\VDA_GMSBMF-main\VDA_25runs\VDA_oMat_AUPR_25runs_rerun.csv');
% 

% writematrix(VDA_AUC, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\VDA_GMSBMF\VDA_GMSBMF-main\VDA_25runs\VDA_LRSSL_AUC_25runs_rerun.csv');
% writematrix(VDA_AUPR, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\VDA_GMSBMF\VDA_GMSBMF-main\VDA_25runs\VDA_LRSSL_AUPR_25runs_rerun.csv');
 
% writematrix(VDA_AUC, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\VDA_GMSBMF\VDA_GMSBMF-main\VDA_25runs\VDA_Cdata_AUC_25runs_rerun.csv');
% writematrix(VDA_AUPR, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\VDA_GMSBMF\VDA_GMSBMF-main\VDA_25runs\VDA_Cdata_AUPR_25runs_rerun.csv');


% writematrix(VDA_AUC, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\VDA_GMSBMF\VDA_GMSBMF-main\VDA_25runs\VDA_Ydata_AUC_5runs_rerun.csv');
% writematrix(VDA_AUPR, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\VDA_GMSBMF\VDA_GMSBMF-main\VDA_25runs\VDA_Ydata_AUPR_5runs_rerun.csv');


writematrix(VDA_AUC, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\VDA_GMSBMF\VDA_GMSBMF-main\VDA_25runs\VDA_hsdn_AUC_25runs_rerun_neworder.csv');
writematrix(VDA_AUPR, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\VDA_GMSBMF\VDA_GMSBMF-main\VDA_25runs\VDA_hsdn_AUPR_25runs_rerun_neworder.csv');
 
% t = toc;

% writematrix(prediction_matrix, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\VDA_GMSBMF\VDA_GMSBMF-main\Figure1run\VDA_hsdn_prediction_rerun.csv');