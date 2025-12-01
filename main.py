import pandas as pd
from sklearn.model_selection import KFold
import os
import sys

sys.path.append("../..")
import numpy as np
import torch
import random
import math
# import argparse
# import json
from model.neural_network import Multi_Omics_Cooperative_Transformer
# import torch.utils.data as Data
from sklearn import preprocessing


def CIndex(pred, ytime_test, ystatus_test):
    N_test = ystatus_test.shape[0]
    ystatus_test = np.squeeze(ystatus_test)
    ytime_test = np.squeeze(ytime_test)
    theta = np.squeeze(pred)
    concord = 0.
    total = 0.
    eav_count = 0
    for i in range(N_test):
        if ystatus_test[i] == 1:
            for j in range(N_test):
                if ytime_test[j] > ytime_test[i]:
                    total = total + 1
                    if theta[j] < theta[i]:
                        concord = concord + 1
                    elif theta[j] == theta[i]:
                        concord = concord + 0.5
                        eav_count = eav_count + 1
                        print("相等的对数")
                        print(eav_count)
    return concord / total


def AUC(pred, ytime_test, ystatus_test):
    N_test = ystatus_test.shape[0]
    ystatus_test = np.squeeze(ystatus_test)
    ytime_test = np.squeeze(ytime_test)
    theta = np.squeeze(pred)
    total = 0
    count = 0

    for i in range(N_test):
        if ystatus_test[i] == 1:
            for j in range(N_test):
                # if ytime_test[i] < quantile_2 < ytime_test[j]:
                if ytime_test[i] < 365 * 1 < ytime_test[j]:
                    total = total + 1
                    if theta[j] < theta[i]:
                        count = count + 1
                    elif theta[j] == theta[i]:
                        count = count + 0.5

    for i in range(N_test):
        if ystatus_test[i] == 1:
            for j in range(N_test):
                # if ytime_test[i] < quantile_2 < ytime_test[j]:
                if ytime_test[i] < 365 * 3 < ytime_test[j]:
                    total = total + 1
                    if theta[j] < theta[i]:
                        count = count + 1
                    elif theta[j] == theta[i]:
                        count = count + 0.5

    for i in range(N_test):
        if ystatus_test[i] == 1:
            for j in range(N_test):
                # if ytime_test[i] < quantile_2 < ytime_test[j]:
                if ytime_test[i] < 365 * 5 < ytime_test[j]:
                    total = total + 1
                    if theta[j] < theta[i]:
                        count = count + 1
                    elif theta[j] == theta[i]:
                        count = count + 0.5

    for i in range(N_test):
        if ystatus_test[i] == 1:
            for j in range(N_test):
                # if ytime_test[i] < quantile_2 < ytime_test[j]:
                if ytime_test[i] < 365 * 10 < ytime_test[j]:
                    total = total + 1
                    if theta[j] < theta[i]:
                        count = count + 1
                    elif theta[j] == theta[i]:
                        count = count + 0.5

    return count / total


if __name__ == '__main__':
    # import time
    # time.sleep(3600*3.5)

    train_rate = 0.8
    train_lr = 2e-4



    base_path = "Cooperative_Transformer_Survival_Analysis_With_Multi_Omics/Datasets/Cancer/"

    standard_scaler = preprocessing.StandardScaler()

    Target_PATH = base_path + "Target/BRCA"
    Target_RNASeq_feature = np.loadtxt(fname=Target_PATH + "/RNASeq_graph_construct.csv", delimiter=",", skiprows=1)[:, 0:6000]
    Target_miRNA_feature = np.loadtxt(fname=Target_PATH + "/miRNA_graph_construct.csv", delimiter=",", skiprows=1)


    Target_ytime = np.loadtxt(fname=Target_PATH + "/ytime.csv", delimiter=",", skiprows=1)
    Target_ystatus = np.loadtxt(fname=Target_PATH + "/ystatus.csv", delimiter=",", skiprows=1)
    Target_RNASeq_feature = standard_scaler.fit_transform(Target_RNASeq_feature)
    Target_miRNA_feature = standard_scaler.fit_transform(Target_miRNA_feature)

    row_num, col_num = Target_miRNA_feature.shape
    for j_col in range(col_num):
        percentiles = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90])
        percentiles_val = np.array([5, 15, 25, 35, 45, 55, 65, 75, 85, 95])
        ptiles_vers = np.percentile(Target_miRNA_feature[:, j_col], percentiles)
        percentiles_val_vers = np.percentile(Target_miRNA_feature[:, j_col], percentiles_val)
        for i_row in range(row_num):

            if Target_miRNA_feature[i_row, j_col] < ptiles_vers[0]:
                Target_miRNA_feature[i_row, j_col] = percentiles_val_vers[0]
            if Target_miRNA_feature[i_row, j_col] >= ptiles_vers[8]:
                Target_miRNA_feature[i_row, j_col] = percentiles_val_vers[9]
            if ptiles_vers[0] <= Target_miRNA_feature[i_row, j_col] < ptiles_vers[1]:
                Target_miRNA_feature[i_row, j_col] = percentiles_val_vers[1]
            if ptiles_vers[1] <= Target_miRNA_feature[i_row, j_col] < ptiles_vers[2]:
                Target_miRNA_feature[i_row, j_col] = percentiles_val_vers[2]
            if ptiles_vers[2] <= Target_miRNA_feature[i_row, j_col] < ptiles_vers[3]:
                Target_miRNA_feature[i_row, j_col] = percentiles_val_vers[3]
            if ptiles_vers[3] <= Target_miRNA_feature[i_row, j_col] < ptiles_vers[4]:
                Target_miRNA_feature[i_row, j_col] = percentiles_val_vers[4]
            if ptiles_vers[4] <= Target_miRNA_feature[i_row, j_col] < ptiles_vers[5]:
                Target_miRNA_feature[i_row, j_col] = percentiles_val_vers[5]
            if ptiles_vers[5] <= Target_miRNA_feature[i_row, j_col] < ptiles_vers[6]:
                Target_miRNA_feature[i_row, j_col] = percentiles_val_vers[6]
            if ptiles_vers[6] <= Target_miRNA_feature[i_row, j_col] < ptiles_vers[7]:
                Target_miRNA_feature[i_row, j_col] = percentiles_val_vers[7]
            if ptiles_vers[7] <= Target_miRNA_feature[i_row, j_col] < ptiles_vers[8]:
                Target_miRNA_feature[i_row, j_col] = percentiles_val_vers[8]

    Target_RNASeq_graph_feature = Target_RNASeq_feature
    Target_miRNA_graph_feature = Target_miRNA_feature

    Target_RNASeq_feature = np.loadtxt(fname=Target_PATH + "/RNASeq_graph_construct.csv", delimiter=",", skiprows=1)[:, 0:6000]
    Target_miRNA_feature = np.loadtxt(fname=Target_PATH + "/miRNA_graph_construct.csv", delimiter=",", skiprows=1)


    Target_ytime = np.loadtxt(fname=Target_PATH + "/ytime.csv", delimiter=",", skiprows=1)
    Target_ystatus = np.loadtxt(fname=Target_PATH + "/ystatus.csv", delimiter=",", skiprows=1)
    Target_RNASeq_feature = standard_scaler.fit_transform(Target_RNASeq_feature)
    Target_miRNA_feature = standard_scaler.fit_transform(Target_miRNA_feature)

    k = 9
    k_assist = 1
    max_cind_list = []
    max_auc_list = []
    for k_num in range(50):
        random.seed(k_num)
        torch.manual_seed(k_num)
        np.random.seed(k_num)
        Target_ystatus_dead_index = np.squeeze(np.argwhere(Target_ystatus == 1))
        Target_ystatus_censor_index = np.squeeze(np.argwhere(Target_ystatus == 0))

        Target_ystatus_dead = Target_ystatus[Target_ystatus_dead_index,]
        Target_ystatus_censor = Target_ystatus[Target_ystatus_censor_index,]

        Target_RNASeq_dead = Target_RNASeq_feature[Target_ystatus_dead_index,]
        Target_RNASeq_censor = Target_RNASeq_feature[Target_ystatus_censor_index,]
        Target_RNASeq_graph_dead = Target_RNASeq_graph_feature[Target_ystatus_dead_index,]
        Target_RNASeq_graph_censor = Target_RNASeq_graph_feature[Target_ystatus_censor_index,]

        Target_miRNA_dead = Target_miRNA_feature[Target_ystatus_dead_index,]
        Target_miRNA_censor = Target_miRNA_feature[Target_ystatus_censor_index,]
        Target_miRNA_graph_dead = Target_miRNA_graph_feature[Target_ystatus_dead_index,]
        Target_miRNA_graph_censor = Target_miRNA_graph_feature[Target_ystatus_censor_index,]

        Target_ytime_dead = Target_ytime[Target_ystatus_dead_index,]
        Target_ytime_censor = Target_ytime[Target_ystatus_censor_index,]

        dead_range = range(Target_RNASeq_dead.shape[0])
        ind_dead_train = random.sample(dead_range, math.floor(Target_RNASeq_dead.shape[0] * train_rate))
        Target_RNASeq_dead_train = Target_RNASeq_dead[ind_dead_train,]
        Target_RNASeq_graph_dead_train = Target_RNASeq_graph_dead[ind_dead_train,]
        print(Target_RNASeq_dead_train.shape)
        Target_miRNA_dead_train = Target_miRNA_dead[ind_dead_train,]
        Target_miRNA_graph_dead_train = Target_miRNA_graph_dead[ind_dead_train,]

        Target_ytime_dead_train = Target_ytime_dead[ind_dead_train,]
        Target_ystatus_dead_train = Target_ystatus_dead[ind_dead_train,]

        ind_dead_rest = [i for i in dead_range if i not in ind_dead_train]
        Target_RNASeq_dead_rest = Target_RNASeq_dead[ind_dead_rest,]
        Target_RNASeq_graph_dead_rest = Target_RNASeq_graph_dead[ind_dead_rest,]

        Target_miRNA_dead_rest = Target_miRNA_dead[ind_dead_rest,]
        Target_miRNA_graph_dead_rest = Target_miRNA_graph_dead[ind_dead_rest,]

        Target_ytime_dead_rest = Target_ytime_dead[ind_dead_rest,]
        Target_ystatus_dead_rest = Target_ystatus_dead[ind_dead_rest,]

        censor_range = range(Target_RNASeq_censor.shape[0])
        ind_censor_train = random.sample(censor_range, math.floor(Target_RNASeq_censor.shape[0] * train_rate))

        Target_RNASeq_censor_train = Target_RNASeq_censor[ind_censor_train,]
        Target_RNASeq_graph_censor_train = Target_RNASeq_graph_censor[ind_censor_train,]

        Target_miRNA_censor_train = Target_miRNA_censor[ind_censor_train,]
        Target_miRNA_graph_censor_train = Target_miRNA_graph_censor[ind_censor_train,]

        Target_ytime_censor_train = Target_ytime_censor[ind_censor_train,]
        Target_ystatus_censor_train = np.squeeze(Target_ystatus_censor[ind_censor_train,])

        Target_RNASeq_train = np.concatenate((Target_RNASeq_dead_train, Target_RNASeq_censor_train), axis=0)
        Target_RNASeq_graph_train = np.concatenate((Target_RNASeq_graph_dead_train, Target_RNASeq_graph_censor_train), axis=0)

        Target_miRNA_train = np.concatenate((Target_miRNA_dead_train, Target_miRNA_censor_train), axis=0)
        Target_miRNA_graph_train = np.concatenate((Target_miRNA_graph_dead_train, Target_miRNA_graph_censor_train), axis=0)
        Target_ystatus_train = np.squeeze(
            np.concatenate((Target_ystatus_dead_train, Target_ystatus_censor_train), axis=0))
        Target_ytime_train = np.squeeze(np.concatenate((Target_ytime_dead_train, Target_ytime_censor_train), axis=0))

        ind_censor_rest = [i for i in censor_range if i not in ind_censor_train]
        Target_RNASeq_censor_rest = Target_RNASeq_censor[ind_censor_rest,]
        Target_RNASeq_graph_censor_rest = Target_RNASeq_graph_censor[ind_censor_rest,]
        Target_miRNA_censor_rest = Target_miRNA_censor[ind_censor_rest,]
        Target_miRNA_graph_censor_rest = Target_miRNA_graph_censor[ind_censor_rest,]
        Target_ytime_censor_rest = np.squeeze(Target_ytime_censor[ind_censor_rest,])
        Target_ystatus_censor_rest = np.squeeze(Target_ystatus_censor[ind_censor_rest,])

        Target_RNASeq_val = np.concatenate((Target_RNASeq_dead_rest, Target_RNASeq_censor_rest), axis=0)
        Target_RNASeq_graph_val = np.concatenate((Target_RNASeq_graph_dead_rest, Target_RNASeq_graph_censor_rest), axis=0)

        Target_miRNA_val = np.concatenate((Target_miRNA_dead_rest, Target_miRNA_censor_rest), axis=0)
        Target_miRNA_graph_val = np.concatenate((Target_miRNA_graph_dead_rest, Target_miRNA_graph_censor_rest), axis=0)

        Target_ytime_val = np.squeeze(np.concatenate((Target_ytime_dead_rest, Target_ytime_censor_rest), axis=0))
        Target_ystatus_val = np.squeeze(np.concatenate((Target_ystatus_dead_rest, Target_ystatus_censor_rest), axis=0))

        Target_RNASeq_val_tensor = torch.tensor(Target_RNASeq_val, dtype=torch.float)
        Target_miRNA_val_tensor = torch.tensor(Target_miRNA_val, dtype=torch.float)
        Target_RNASeq_train_tensor = torch.tensor(Target_RNASeq_train, dtype=torch.float)
        Target_miRNA_train_tensor = torch.tensor(Target_miRNA_train, dtype=torch.float)

        Target_RNASeq_graph_val_tensor = torch.tensor(Target_RNASeq_graph_val, dtype=torch.float)
        Target_miRNA_graph_val_tensor = torch.tensor(Target_miRNA_graph_val, dtype=torch.float)
        Target_RNASeq_graph_train_tensor = torch.tensor(Target_RNASeq_graph_train, dtype=torch.float)
        Target_miRNA_graph_train_tensor = torch.tensor(Target_miRNA_graph_train, dtype=torch.float)

        Target_RNASeq_graph_tensor = torch.cat((Target_RNASeq_graph_train_tensor, Target_RNASeq_graph_val_tensor), axis=0)
        Target_RNASeq_graph_tensor1 = torch.unsqueeze(Target_RNASeq_graph_tensor, 1)  # N*1*d
        Target_RNASeq_graph_tensor2 = torch.unsqueeze(Target_RNASeq_graph_tensor, 0)  # 1*N*d
        Target_W_Gene_distance = ((Target_RNASeq_graph_tensor1 - Target_RNASeq_graph_tensor2) ** 2).sum(2)  # N*N*d -> N*N
        Target_W_Gene_distance_temp = Target_W_Gene_distance.reshape(-1, 1)
        Target_distance = torch.median(Target_W_Gene_distance_temp, 0)
        Target_W_Gene_cof = torch.exp(-Target_W_Gene_distance / (0.3 *Target_distance[0]))
        Target_W_Gene = torch.zeros_like(Target_W_Gene_cof)
        Target_W_Gene_assist = torch.zeros_like(Target_W_Gene_cof)
        if k > 0:
            topk, indices = torch.topk(Target_W_Gene_cof, k)
            # print(indices)
            mask = torch.zeros_like(Target_W_Gene_cof)
            mask = mask.scatter(1, indices, 1)
            # mask = ((mask + torch.t(mask)) > 0).type(torch.float32)  # union, kNN graph
            # mask = ((mask > 0) & (torch.t(mask) > 0)).type(torch.float32)  # intersection, kNN graph
            Target_W_Gene = mask

            topk_assist, indices_assist = torch.topk(Target_W_Gene_cof, k_assist)
            # print(indices)
            mask = torch.zeros_like(Target_W_Gene_cof)
            mask = mask.scatter(1, indices_assist, 1)
            # mask = ((mask + torch.t(mask)) > 0).type(torch.float32)  # union, kNN graph
            # mask = ((mask > 0) & (torch.t(mask) > 0)).type(torch.float32)  # intersection, kNN graph
            Target_W_Gene_assist = mask
        # W_Gene = W_Gene / W_Gene.sum(0)
        # print(W_Gene)
        Target_miRNA_graph_tensor = torch.cat((Target_miRNA_graph_train_tensor, Target_miRNA_graph_val_tensor), axis=0)
        Target_miRNA_graph_tensor1 = torch.unsqueeze(Target_miRNA_graph_tensor, 1)  # N*1*d
        Target_miRNA_graph_tensor2 = torch.unsqueeze(Target_miRNA_graph_tensor, 0)  # 1*N*d
        Target_W_miRNA_distance = ((Target_miRNA_graph_tensor1 - Target_miRNA_graph_tensor2) ** 2).sum(2)  # N*N*d -> N*N
        Target_W_miRNA_distance_temp = Target_W_miRNA_distance.reshape(-1, 1)
        Target_distance = torch.median(Target_W_miRNA_distance_temp, 0)
        Target_W_miRNA_cof = torch.exp(-Target_W_miRNA_distance / (0.2 *Target_distance[0]))
        Target_W_miRNA = torch.zeros_like(Target_W_miRNA_cof)
        Target_W_miRNA_assist = torch.zeros_like(Target_W_miRNA_cof)
        if k > 0:
            topk, indices = torch.topk(Target_W_miRNA_cof, k)
            mask = torch.zeros_like(Target_W_miRNA_cof)
            mask = mask.scatter(1, indices, 1)
            # mask = ((mask + torch.t(mask)) > 0).type(torch.float32)  # union, kNN graph
            # mask = ((mask > 0) & (torch.t(mask) > 0)).type(torch.float32)  # intersection, kNN graph
            Target_W_miRNA = mask

            topk_assist, indices_assist = torch.topk(Target_W_miRNA_cof, k_assist)
            # print(indices)
            mask = torch.zeros_like(Target_W_miRNA_cof)
            mask = mask.scatter(1, indices_assist, 1)
            # mask = ((mask + torch.t(mask)) > 0).type(torch.float32)  # union, kNN graph
            # mask = ((mask > 0) & (torch.t(mask) > 0)).type(torch.float32)  # intersection, kNN graph
            Target_W_miRNA_assist = mask

        # W_miRNA = W_miRNA / W_miRNA .sum(0)
        # print(W_miRNA)
        # Target_W = (Target_W_Gene_cof + Target_W_miRNA_cof)
        eps = np.finfo(float).eps
        W = (Target_W_Gene_cof + Target_W_miRNA_cof) / 2 + torch.eye(Target_W_Gene_cof.shape[0])
        D = torch.sum(W, dim=1)
        D_sqrt_inv = torch.sqrt(1.0 / (D + eps))
        Target_W = D_sqrt_inv * W * D_sqrt_inv

        Target_W_Temp = ((Target_W_Gene > 0) | (Target_W_miRNA > 0)).type(torch.float32)  # intersection, kNN graph
        Target_W_assist = ((Target_W_Gene_assist > 0) & (Target_W_miRNA_assist > 0)).type(torch.float32)  # intersection, kNN graph
        Target_W_Trans = torch.matmul(Target_W_Temp , Target_W_assist)
        Target_W_Trans[Target_W_Trans > 0] = 1

        model = Multi_Omics_Cooperative_Transformer(p_RNASeq=0.5, p_miRNA=0.5, p_graph=0)

        optimizer = torch.optim.Adam([{'params': model.parameters()}, ], lr=train_lr, weight_decay=5e-4)
        max_cind = 0.0
        max_auc = 0.0
        best_iter = 0
        for iter in range(100):
            model.train()

            index = np.squeeze(np.arange(0, Target_RNASeq_train.shape[0]))
            fir_index = np.random.choice(index, size=Target_RNASeq_train.shape[0], replace=False)
            # print(fir_index)
            # seco_index = np.array(list(set(index.tolist()) - set(fir_index)))

            Target_ystatus_batch_train = Target_ystatus_train[fir_index,]
            Target_ystatus_train_tensor = torch.tensor(Target_ystatus_batch_train, dtype=torch.float)

            Target_ytime_batch_train = Target_ytime_train[fir_index,]
            Target_ytime_train_tensor = torch.tensor(Target_ytime_batch_train, dtype=torch.float)

            Target_RNASeq_tensor = torch.cat((Target_RNASeq_train_tensor, Target_RNASeq_val_tensor), 0)
            Target_miRNA_tensor = torch.cat((Target_miRNA_train_tensor, Target_miRNA_val_tensor), 0)
            # DNA_tensor = torch.cat((DNA_train_tensor, DNA_val_tensor), 0)

            theta = model.get_survival_result(Target_RNASeq_tensor, Target_miRNA_tensor,
                                                     Target_W, Target_W_Trans)[fir_index,]

            real_batch_size = Target_ystatus_train_tensor.shape[0]
            R_matrix_batch_train = torch.tensor(np.zeros([real_batch_size, real_batch_size], dtype=int),
                                                dtype=torch.float)

            for i in range(real_batch_size):
                R_matrix_batch_train[i,] = torch.tensor(
                    np.array(list(map(int, (Target_ytime_train_tensor >= Target_ytime_train_tensor[i])))))
            model.train()
            exp_theta = torch.reshape(torch.exp(theta), [real_batch_size])
            theta = torch.reshape(theta, [real_batch_size])
            fuse_loss = -torch.mean(torch.mul((theta - torch.log(torch.sum(torch.mul(exp_theta,
                                                                                     R_matrix_batch_train),
                                                                           dim=1))),
                                              torch.reshape(Target_ystatus_train_tensor, [real_batch_size])))

            loss = fuse_loss

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            eval_model = Multi_Omics_Cooperative_Transformer(p_RNASeq=0, p_miRNA=0, p_graph=0)
            eval_model.load_state_dict(model.state_dict())  # copy? looks okay
            eval_model.eval()
            pred_train = eval_model.get_survival_result(Target_RNASeq_tensor, Target_miRNA_tensor,
                                                        Target_W, Target_W_Trans)[0:Target_RNASeq_train.shape[0], ]

            cind_train = CIndex(pred_train, Target_ytime_train, Target_ystatus_train)
            # auc_train = AUC(theta, ytime_train, ystatus_train)

            with open('Cox_Result_lr5e-05/Cox_Train/Regularizer_eval_randomseed' + str(k_num) + '.log',
                      'a') as f:
                f.writelines('Iteration:' + str(iter) + ",cind_train:" + str(cind_train) + '\n')
            pred_val = eval_model.get_survival_result(Target_RNASeq_tensor, Target_miRNA_tensor,
                                                     Target_W, Target_W_Trans)[Target_RNASeq_train.shape[0]:, :]
            cind_val = CIndex(pred_val.detach().numpy(), Target_ytime_val, Target_ystatus_val)
            auc_val = AUC(pred_val.detach().numpy(), Target_ytime_val, Target_ystatus_val)

            with open('Cox_Result_lr5e-05/Cox_Val/Regularizer_eval_randomseed' + str(k_num) + '.log',
                      'a') as f:
                f.writelines(
                    'Iteration:' + str(iter) + ",cind:" + str(cind_val) + ",auc:" + str(auc_val) + '\n')
            if cind_train - cind_val > 0.05:
                if cind_val >= max_cind:
                    max_cind = cind_val
                    max_auc = auc_val
                    best_iter = iter
                    count = 0
                else:
                    count = count + 1
            if iter > 100:
                break

        max_cind_list.append(max_cind)
        max_auc_list.append(max_auc)
        with open('Cox_Result_lr5e-05/Cox_Val/Regularizer_eval_AVE_Cind' + '.log', 'a') as f:
            f.writelines(str(best_iter) + ',' + str(max_cind) + ',' + str(max_auc) + '\n')
    with open('Cox_Result_lr5e-05/Cox_Val/Regularizer_eval_AVE_Cind' + '.log', 'a') as f:
        f.writelines("max cind ave:" + str(np.mean(np.array(max_cind_list))) + '\n')
        f.writelines("max cind std:" + str(np.std(np.array(max_cind_list))) + '\n')
        f.writelines("max auc ave:" + str(np.mean(np.array(max_auc_list))) + '\n')
        f.writelines("max auc std:" + str(np.std(np.array(max_auc_list))) + '\n')
