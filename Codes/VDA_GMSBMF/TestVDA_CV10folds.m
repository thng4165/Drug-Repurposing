clear all
rng('default'); seed = 12345; rng(seed)
tic;
% DirOut = ['Results']; if ~exist(DirOut,'dir'); mkdir(DirOut); end 
% test = '' ;

%datadirSet = {'VDdataset1'};
methodset  = {'VDA-GMSBMF'};
%for i_datadir = 1:length(datadirSet)
    %datadir = datadirSet{i_datadir}; 
    nCV = 5; nfold  = 10; CVtype = 'CVa';    
  
    %% 1. Load Datesets
     load Datasets\MatlabDataFiles\HDVD    
    % load Datasets\MatlabDataFiles\LAGCN
    % load Datasets\MatlabDataFiles\Fdataset
    % load Datasets\MatlabDataFiles\Cdataset
    % load Datasets\MatlabDataFiles\LRSSL
    % load Datasets\MatlabDataFiles\Ydataset
   

     Wdv = didr'; 
     Wvv = disease;
     Wdd = drug; 

    % oMat-MechDB dataset 
    % a = load('Datasets\MatlabDataFiles\rare_disease_drug.mat');
    % Wdv = a.data;   
    % b = load('Datasets\MatlabDataFiles\rare_drug_sim.mat');
    % Wdd = b.data;
    % c = load('Datasets\MatlabDataFiles\rare_disease_sim.mat');
    % Wvv = c.data;

    % hsdn-MechDB dataset
    % a = load('Datasets\MatlabDataFiles\hsdn_MechDB_dd_association_numeric.mat');
    % Wdv = a.data; 
    % b = load('Datasets\MatlabDataFiles\hsdn_MechDB_drug_sim.mat');
    % Wdd = b.data;
    % c = load('Datasets\MatlabDataFiles\hsdn_MechDB_disease_sim_GIP.mat');
    % Wvv = c.data;



    [Nd,Nv] = size(Wdv); 
    % % %  AUROC_set = [];AUPRC_set =[]; Acc_set = []; Sen_set = []; Spe_set = [];Pre_set =[];  
    % tbScalarMeanStdAll = table;  
    % tbVecMeanStdAll    = table; 
    % 
	% tbScalarMeanAll = table ;
	% tbScalarStdAll  = table ;
	% tbVecMeanAll    = table ;
	% tbVecStdAll     = table ;
 
    for i_mtd = 1:length( methodset )
        method      = methodset{i_mtd} ;
        MatVecSet   = []; 
        tbScalarSet = [];   
        ii_fold = 0; 
		%disp([datadir,'--',method,': i_cv-from 1 to ',num2str(nCV),'; i_fold-from 1 to ',num2str(nfold)])
        %for i_cv = 1:nCV 
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
                                gm= 0.5; w=0.4;lambda1 = 0.5; lambda2 = lambda1; lambda3 = lambda2;  tol1 = 2*1e-30;tol2 = 1*1e-40; 
                        %end                          
                        [M_recovery] = A_VDA_GMSBMF(matDV, Sdd, Svv, gm,w, lambda1, lambda2, lambda3, floor( min(size(matDV)) * 1), tol1, tol2, 400) ; 
                            
                    otherwise; error(['There is no defintion: ', method]);                  
                end 

                % 
                MatPredict(Ind_test)= M_recovery(Ind_test);   
                % labels = Wdv(Ind_test) ; 
                % scores = M_recovery(Ind_test); 
                % %
                % [tbScalar, tbVec]=getPerfMetricROCcompute( scores,labels,1, 0);   
                % ii_fold = ii_fold + 1 ; 
                % tbScalar.Properties.RowNames = strcat([method,'_', num2str(ii_fold),'_'],tbScalar.Properties.RowNames ); 
                % tbScalarSet{ii_fold}  = tbScalar ;
                % StrMetricVec          = strcat( [method,'_'], tbVec.Properties.RowNames  );  
                % MatVecSet{ii_fold}    = tbVec{:,:} ;
                % 
            end 

            % matRealLabel = Wdv; 
            % np = nnz(matRealLabel);
            % vv = MatPredict(:); 
            % [vs,ind] = sort(vv, 'descend' );   
            % matPredLabel = false( size(matRealLabel) )  ;  
            % matPredLabel( ind(1:np)  ) = 1;   
            % % matRealLabel     matPredLabel        

        %end
        % %   
        % tbScalarSet1 = cat(1, tbScalarSet{:});  
        % tbScalarMean = array2table( mean(tbScalarSet1{:,:},1), 'RowNames',{[method,'_mean'  ]}, 'VariableNames', tbScalarSet1.Properties.VariableNames );  
        % tbScalarStd  = array2table( std(tbScalarSet1{:,:},1), 'RowNames',{[method,'_std'  ]}  , 'VariableNames', tbScalarSet1.Properties.VariableNames );  
        % % tbScalarMeanAll = [tbScalarMeanAll; tbScalarMean] ;
        % % tbScalarStdAll  = [tbScalarStdAll; tbScalarStd] ;
        % tbScalarMeanAll = [ tbScalarMeanAll; tbScalarMean ] ; 
        % tbScalarStdAll  = [ tbScalarStdAll; tbScalarStd ]  ;
        %
        % RowVecSet       = cat(3, MatVecSet{:}); % RowVec = mean(RowVec,3);  
        % tbVecMean    = array2table(mean(RowVecSet,3), 'RowNames',  strcat(['mean_'],StrMetricVec )  );  
        % tbVecStd     = array2table(std(RowVecSet,[], 3), 'RowNames',  strcat(['std_'],StrMetricVec )   ); 
        % tbVecMeanAll = [tbVecMeanAll; [tbVecMean   ]];  
        % tbVecStdAll = [tbVecStdAll;   [tbVecStd    ]];  

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

time = toc;
