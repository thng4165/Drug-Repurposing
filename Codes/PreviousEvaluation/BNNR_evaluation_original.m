
   
myseed = 2024;
rng(myseed);
    %% 1. Load Datasets%%%%%%%%%%%%%%%%%%%%%%%%
    % Load the dataset specified in the input
    % load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\Fdataset
    % load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\Cdataset
    load C:\Trang\KIProjects\ComprehensionDR\Datasets\MatlabDataFiles\Ydataset
    
    % Extract matrices from loaded data
    Wrr = drug;
    Wdd = disease;
    Wdr = didr;
    Wrd = Wdr';
    A = Wrd;

    folds = 10;
    
    % Split the data into training and testing sets
    [positiveId, crossval_id] = train_test_split(A, folds, '1');
    O = [];
    P = [];
  for fold = 1:folds
        train = A;
        PtestID = positiveId(crossval_id == fold);
        negativeID = find(Wrd == 0);
        num = numel(negativeID);
        Nidx = randperm(num);
        NtestID = negativeID(Nidx(1:length(PtestID)));
        train(PtestID) = 0;
        Wdr_train = train';
        [dn, dr] = size(Wdr_train);


        %%%%%%%% 2. Algorithm Code %%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Algorithm Parameters
        maxiter = 300;
        alpha = 1;
        beta = 10;
        tol1 = 2 * 1e-3;
        tol2 = 1 * 1e-5;
        
        % Construct the combined matrix T
        T = [Wrr, Wdr_train'; Wdr_train, Wdd];
        [t1, t2] = size(T);
        
        % Create a binary mask for T
        trIndex = double(T ~= 0);
        
        % Apply the BNNR algorithm to recover missing values in T
        [WW, iter] = BNNR(alpha, beta, T, trIndex, tol1, tol2, maxiter, 0, 1);
        
        % Extract the recovered matrix corresponding to A
        M_recovery = WW((t1 - dn + 1) : t1, 1 : dr);
        Rt = M_recovery';

        % Prepare data for evaluation
        origin = [A(PtestID); A(NtestID)];
        pred = [Rt(PtestID); Rt(NtestID)];
        O = [O;origin];
        P = [P;pred];
   
end
 inputObs_matrix = O; % transpose to sort by column - disease
 prediction_matrix = P;



[roc_x, roc_y, ~, AUC] = perfcurve(inputObs_matrix, prediction_matrix, true);

[pr_x, pr_y, ~, AUPR] = perfcurve(inputObs_matrix, prediction_matrix, 1, 'xCrit', 'reca', 'yCrit', 'prec');

writematrix(roc_x, 'roc_x_Ydata.csv');
writematrix(roc_y, 'roc_y_Ydata.csv');
writematrix(pr_x, 'pr_x_Ydata.csv');
writematrix(pr_y, 'pr_y_Ydata.csv');

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