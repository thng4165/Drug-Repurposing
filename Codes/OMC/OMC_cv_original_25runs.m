clear all
addpath('C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\OMC_method\OMC-master');

%% 1. Load Datesets

    % load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\Fdataset
     % load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\Cdataset
    load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\Ydataset
    

Wrr = drug;
Wdd = disease;
Wdr = didr;
Wrd = Wdr';
A = Wdr;

OMC_AUC = [];
OMC_AUPR = [];
seed_values = randsample(1000000:2000000, 16, false); 
parfor iseed = 1:length(seed_values)
 
    fprintf('iseed: %d\n', iseed);

%% Cross Validation


    nfolds = 10;
    myseed = seed_values(iseed);
    rng(myseed, 'twister');

    [positiveId, crossval_id] = train_test_split(A, nfolds, '1');
    O = [];
    P = [];

   
%% 2. OMC2 algorithm
alpha = 1;
beta = 10;
K = 10;
tol1 = 2*1e-3;
tol2 = 1*1e-5;
maxiter = 300;
[nr,nc] = size(Wdr);

 
for fold = 1:nfolds
   fold
        PtestID = positiveId(crossval_id == fold);
        negativeID = find(A == 0);
        num = numel(negativeID);
        Nidx = randperm(num);
        NtestID = negativeID(Nidx(1:length(PtestID)));
           
    test_idx  = PtestID;

    % train_data = A; 
    % train_data(test_idx) = NaN;
    
    P_TMat = A;
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

    % Prepare data for evaluation
        origin = [A(PtestID); A(NtestID)];
        pred = [M_recovery(PtestID); M_recovery(NtestID)];
        O = [O;origin];
        P = [P;pred];
end

%%  CV: AUC-AUPR
 inputObs_matrix = O; % transpose to sort by column - disease
 prediction_matrix = P;



[roc_x, roc_y, ~, AUC] = perfcurve(inputObs_matrix, prediction_matrix, true);

[pr_x, pr_y, ~, AUPR] = perfcurve(inputObs_matrix, prediction_matrix, 1, 'xCrit', 'reca', 'yCrit', 'prec');

OMC_AUC(iseed) = AUC;
 OMC_AUPR(iseed) = AUPR;
end

writematrix(OMC_AUC, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\OMC_method\OMC-master\Pr_roc_CV_original\OMC_original_Ydata_AUC_16runs_rerun.csv');
writematrix(OMC_AUPR, 'C:\Trang\KIProjects\ComprehensionDR\RerunMatlabCode\OMC_method\OMC-master\Pr_roc_CV_original\OMC_original_Ydata_AUPR_16runs_rerun.csv');

% writematrix(roc_x, 'OMC_roc_x_Ydata.csv');
% writematrix(roc_y, 'OMC_roc_y_Ydata.csv');
% writematrix(pr_x, 'OMC_pr_x_Ydata.csv');
% writematrix(pr_y, 'OMC_pr_y_Ydata.csv');
% 
% % Plot ROC curve
% figure;
% plot(roc_x, roc_y, 'b-', 'LineWidth', 2);
% xlabel('False Positive Rate','FontSize', 12);
% ylabel('True Positive Rate','FontSize', 12);
% title(['ROC Curve (AUC = ' num2str(AUC) ')'], 'FontSize', 14);
% grid on;
% hold on;
% plot([0, 1], [0, 1], '--', 'Color', [0.5, 0.5, 0.5]);
% hold off;
% 
% 
% figure;
% plot(pr_x, pr_y, 'b-', 'LineWidth', 2);
% xlabel('Recall', 'FontSize', 12);
% ylabel('Precision', 'FontSize', 12);
% title(['Precision-Recall Curve (AUPR = ' num2str(AUPR) ')'], 'FontSize', 14);
% grid on;
% 
