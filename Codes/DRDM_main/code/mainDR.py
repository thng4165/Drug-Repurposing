#import time
#start_time = time.time()
import argparse
from loaderDR import *
import torch
from torch.optim.lr_scheduler import CyclicLR
from model import DRDM
from sklearn import metrics
import numpy as np
from collections import OrderedDict

import random
import sys

sys.path.append(".")
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

def parse_args():
    parser = argparse.ArgumentParser(description="Run DRDM.")
    parser.add_argument('--dataset', nargs='?', default='Fdataset', help='Choose a dataset. [Fdataset/Cdataset/LRSSL/HDVD/LAGCN/Ydataset/oMatMechDB/hsdnMechDB/Fdataset_mul/Cdataset_mul/Ydataset_mul]')
    parser.add_argument('--epochs', type=int, default=200, help='Number of epochs.')
    parser.add_argument('--batch_size', type=int, default=1024*5, help='Batch size.')
    parser.add_argument('--lr', type=float, default=0.01, help='init Learning rate.')
    parser.add_argument('--embedding_size', type=int, default=64)
    parser.add_argument('--n_layers', type=int, default=2)
    parser.add_argument('--weight_decay', type=float, default=0.0001)
    parser.add_argument('--disease_TopK', type=int, default=4)
    parser.add_argument('--drug_TopK', type=int, default=4)
    parser.add_argument("--seed", type=int, default=666)
    parser.add_argument("--n_splits", type=int, default=10)
    parser.add_argument("--num_trials", type=int, default=9)
    parser.add_argument('--exp_coff', type=float, default=0.9)
    parser.add_argument('--wd', type=float, default=0.6, help='the coefficient of feature fusion ')
    parser.add_argument('--wr', type=float, default=0.6, help='the coefficient of feature fusion ')
    return parser.parse_args()

def config_model():
    config = OrderedDict()
    config['dataset'] = args.dataset
    config['epochs'] = args.epochs
    config['batch_size'] = args.batch_size
    config['lr'] = args.lr
    config['embedding_size'] = args.embedding_size
    config['n_layers'] = args.n_layers
    config['weight_decay'] = args.weight_decay
    config["disease_TopK"] = args.disease_TopK
    config['drug_TopK'] = args.drug_TopK
    config['seed'] = args.seed
    config['n_splits'] = args.n_splits
    config['device'] = device
    config['exp_coff'] = args.exp_coff
    config['wd'] = args.wd
    config['wr'] = args.wr
    return config

def setup_seed(seed):
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    np.random.seed(seed)
    random.seed(seed)

def sort_matrix(score_matrix, label_matrix):
    score_matrix = score_matrix.T
    label_matrix = label_matrix.T
    # For each column, sort the scores and reorder the labels accordingly
    sort_index = np.zeros_like(score_matrix, dtype=int)
    score_sorted = np.zeros_like(score_matrix)
    label_sorted = np.zeros_like(label_matrix)
    for col in range(score_matrix.shape[1]):
        # Sort descending (highest score first)
        order = np.argsort(-score_matrix[:, col])
        sort_index[:, col] = order
        score_sorted[:, col] = score_matrix[order, col]
        label_sorted[:, col] = label_matrix[order, col]
    return score_sorted, label_sorted, sort_index

def compute_curve_metrics(label_sorted):
    # Compute TPR, FPR, Recall, Precision for each cutoff (row)
    n_rows, n_cols = label_sorted.shape
    tpr_list, fpr_list, recall_list, precision_list = [], [], [], []
    for cutoff in range(1, n_rows + 1):
        P_matrix = label_sorted[:cutoff, :]
        N_matrix = label_sorted[cutoff:, :]
        TP = np.sum(P_matrix > 0)
        FP = np.sum(P_matrix == 0)
        TN = np.sum(N_matrix == 0)
        FN = np.sum(N_matrix > 0)
        tpr = TP / (TP + FN) if (TP + FN) > 0 else 0
        fpr = FP / (FP + TN) if (FP + TN) > 0 else 0
        recall = tpr
        precision = TP / (TP + FP) if (TP + FP) > 0 else 0
        tpr_list.append(tpr)
        fpr_list.append(fpr)
        recall_list.append(recall)
        precision_list.append(precision)
    return np.array(tpr_list), np.array(fpr_list), np.array(recall_list), np.array(precision_list)

def auc_aupr_from_curve(fpr_list, tpr_list, recall_list, precision_list):
    # Use numpy.trapz for trapezoidal integration
    auc = np.trapz(tpr_list, fpr_list)
    aupr = np.trapz(precision_list, recall_list)
    return auc, aupr

args = parse_args()
config = config_model()

if __name__ == "__main__":
    avg_auroc, avg_aupr = [], []
    avg_auc_curve, avg_aupr_curve = [], []
    for i in range(args.num_trials):
        setup_seed(i+10)
        disease_adj, drug_adj, original_interactions, all_train_mask, all_test_mask, pos_weight = data_preparation(args)
        all_scores, all_labels, all_indices = [], [], []
        print(f'+++++++++++++++This is {i + 1}-th 10 fold validation.+++++++++++++++')
        for fold_num in range(len(all_train_mask)):
            print(f'---------------This is {fold_num + 1}-th fold validation.---------------')
            # dataset splitting
            train_manager, test_manager = data_split(config, all_train_mask[fold_num], all_test_mask[fold_num],
                                                     original_interactions)
            train_adj = train_manager.train_adj
            # model loading and initialization
            model = DRDM(config, (train_manager, train_adj, disease_adj, drug_adj, pos_weight)).to(device)
            # training
            optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
            lr_scheduler = CyclicLR(optimizer,
                                    base_lr=0.1 * args.lr,
                                    max_lr=args.lr,
                                    step_size_up=20,
                                    mode="exp_range",
                                    gamma=0.995,
                                    cycle_momentum=False)
            for epoch in range(args.epochs):
                model.train()
                loss_list = []
                for batch in train_manager.iter_batch(shuffle=True):
                    diseases, drugs, labels, _ = batch  # Only use the first three for training
                    loss, scores = model.forward([diseases, drugs, labels], True)
                    model.zero_grad()
                    loss.backward()
                    optimizer.step()
                    lr_scheduler.step()
                model.eval()
                scores, labels, indices = [], [], []
                for batch in test_manager.iter_batch():
                    # Unpack batch to get indices
                    if len(batch) == 4:
                        input_disease, input_drug, label, batch_indices = batch
                    else:
                        input_disease, input_drug, label = batch
                        batch_indices = list(zip(input_disease, input_drug))
                    score, label = model.predict([input_disease, input_drug, label])
                    scores.append(score.cpu().detach().numpy())
                    labels.append(label)
                    indices.extend(batch_indices)
                loss_sum = np.sum(loss_list)
                scores = np.concatenate(scores)
                labels = np.concatenate(labels)
                # aupr = metrics.average_precision_score(y_true=labels, y_score=scores)
                # auroc = metrics.roc_auc_score(y_true=labels, y_score=scores)
                # print(f'Epoch: {epoch + 1}, auroc: {auroc}, aupr: {aupr}')
                if (epoch + 1) == args.epochs:
                    all_scores.append(scores)
                    all_labels.append(labels)
                    all_indices.extend(indices)
                    # Print predicted scores and labels (first 10 values)
                    # print(f">> Fold {fold_num + 1} predicted scores (first 10): {scores[:10]}")
                    # print(f">> Fold {fold_num + 1} true labels (first 10): {labels[:10]}")

                    # Save predictions and labels for this fold
                    # np.save(f"fold_{i + 1}_fold_{fold_num + 1}_scores.npy", scores)
                    # np.save(f"fold_{i + 1}_fold_{fold_num + 1}_labels.npy", labels)
                    # np.save(f"fold_{i + 1}_fold_{fold_num + 1}_indices.npy", np.array(indices))
        all_scores = np.concatenate(all_scores)
        all_labels = np.concatenate(all_labels)
        
        # Reconstruct matrices with original positions
        all_scores_matrix = np.zeros_like(original_interactions, dtype=np.float32)
        all_labels_matrix = np.zeros_like(original_interactions, dtype=np.float32)
        for idx, (i_pos, j_pos) in enumerate(all_indices):
            all_scores_matrix[i_pos, j_pos] = all_scores[idx]
            all_labels_matrix[i_pos, j_pos] = all_labels[idx]

        # --- Column-wise sorting and curve calculation ---
        score_sorted, label_sorted, sort_index = sort_matrix(all_scores_matrix, all_labels_matrix)
       

        tpr_list, fpr_list, recall_list, precision_list = compute_curve_metrics(label_sorted)
        auc_curve, aupr_curve = auc_aupr_from_curve(fpr_list, tpr_list, recall_list, precision_list)
        avg_auc_curve.append(auc_curve)
        avg_aupr_curve.append(aupr_curve)
        print(f"Column-sorted curve AUC: {auc_curve:.5f}, AUPR: {aupr_curve:.5f}")

        print(f"Column-sorted curve AUC: {avg_auc_curve}, AUPR: {avg_aupr_curve}")