from itertools import *
import math
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn import preprocessing


class Multi_Omics_Cooperative_Transformer(nn.Module):

    def __init__(self, p_RNASeq=0, p_miRNA=0, p_graph=0):
        nn.Module.__init__(self)

        self.RNASeq_Wq = nn.Sequential(nn.Linear(6000, 100), nn.BatchNorm1d(100), nn.Tanh())

        self.miRNA_Wq = nn.Sequential(nn.Linear(600, 100), nn.BatchNorm1d(100), nn.Tanh())

        self.dna_Wq = nn.Sequential(nn.Linear(6000, 100), nn.BatchNorm1d(100), nn.Tanh())

        self.Wq_1 = nn.Linear(100, 50)
        self.Wk_1 = nn.Linear(100, 50)
        self.Wv_1 = nn.Linear(100, 50)

        self.Wq_2 = nn.Linear(100, 50)
        self.Wk_2 = nn.Linear(100, 50)
        self.Wv_2 = nn.Linear(100, 50)

        self.atten_fuse_weight = nn.Linear(100, 50)
        self.BN = nn.BatchNorm1d(100)
        self.Graph_BN = nn.BatchNorm1d(100)
        self.p_RNASeq = p_RNASeq
        self.p_miRNA = p_miRNA
        self.p_graph = p_graph
        self.Graph_Wq_value = nn.Linear(100, 100)
        self.Graph_Wq_query = nn.Linear(100, 100)
        self.Graph_Wq_key = nn.Linear(100, 100)

        self.bm_full = nn.Sequential(nn.Linear(100, 30, bias=False), nn.ReLU())
        self.Cox_layer = nn.Linear(30, 1, bias=False)


    def get_gene_feature(self, RNASeq_feature):
        RNASeq_feature = self.RNASeq_Wq(RNASeq_feature)
        return RNASeq_feature

    def get_miRNA_feature(self, miRNA_feature):
        miRNA_feature = self.miRNA_Wq(miRNA_feature)
        return miRNA_feature

    def get_graph_feature(self, fused_feature, S, S_Trans):
        phi_fea = torch.matmul(S, fused_feature) / torch.reshape(torch.sum(S, dim=1), (-1, 1))
        # Target_S_Temp = (S > 0).type(torch.float32)
        # phi_fea = torch.matmul(Target_S_Temp, fused_feature) / torch.reshape(torch.sum(Target_S_Temp, dim=1), (-1, 1))
        fused_feature = F.dropout(fused_feature, self.p_graph)
        value_1 = self.Graph_Wq_value(fused_feature).unsqueeze(0)
        query_1 = (self.Graph_Wq_query(fused_feature)).unsqueeze(0)
        key_1 = (self.Graph_Wq_key(fused_feature)).unsqueeze(0)

        # 求注意力分数
        attention_weights = torch.matmul(query_1, query_1.transpose(-2, -1)) / (100 ** 0.5)
        # attention_weights = torch.matmul(query_1, key_1.transpose(-2, -1)) / (100 ** 0.5)
        # 对注意力分数归一化，最后一个维度做
        # if torch.isnan(attention_weights).any():
        #     print("W_graph_1_temp contains NaN values.")
        # if torch.isinf(attention_weights).any() or torch.isinf(torch.exp(attention_weights)).any():
        #     print("One of the input tensors contains infinite values.")
        # print(attention_weights)
        if torch.isinf(attention_weights).any():
            print("One of the input tensors contains infinite values.")
        if torch.isinf(torch.exp(attention_weights)).any():
            print("One of the input tensors contains infinite values.")
        print(attention_weights.shape)
        print(S_Trans.shape)
        S_norm = S 
        attention_weights_1 = torch.exp(attention_weights ) * S_Trans * S_norm/ torch.sum(torch.exp(attention_weights)* S_norm *
                                                                                 S_Trans, dim=-1).reshape(-1, 1)
        # attention_weights = nn.functional.softmax(attention_weights, dim=-1)
        # 将归一化后的权重与value相乘

        # attention_weights_2 = (attention_weights * S_Trans * S_norm) / torch.sum(attention_weights * S_Trans * S_norm, dim=-1).reshape(-1, 1)
        # attention_weights_2 = (attention_weights_1 * S_norm) / torch.sum(attention_weights_1 * S_norm,
        #                                                                  dim=-1).reshape(-1, 1)
        # attention_weights_2 = (attention_weights_1 + S_Trans * S_norm) 
        attention_weights_2 = attention_weights_1 
        attended_values = torch.matmul(torch.squeeze(attention_weights_2),
                                       F.normalize(torch.squeeze(value_1)))
        # attended_values = torch.cat((torch.squeeze(torch.matmul(S_norm, value_1)), self.Graph_Wq_value(fused_feature)), 1)
        # attended_values = torch.cat((F.normalize(torch.squeeze(torch.matmul(attention_weights_2, value_1))), F.normalize(self.Graph_Wq_value(fused_feature))), 1)
        # attended_values = torch.cat((torch.squeeze(F.normalize(torch.matmul(S_norm, value_1))), F.normalize(self.Graph_Wq_value(fused_feature))), 1)
        # attended_values = torch.squeeze(torch.matmul(S, value_1))
        return self.bm_full(self.Graph_BN(attended_values))

    def get_Transformer_feature(self, gene, miRNA):
        # x: [batch_size, seq_len, input_dim]
        gene = torch.squeeze(gene)
        miRNA = torch.squeeze(miRNA)

        gene_1 = gene.unsqueeze(1)
        miRNA_1 = miRNA.unsqueeze(1)
        x = torch.cat((gene_1, miRNA_1), 1)
        print(x.shape)
        batch_size, seq_len, _ = x.size()
        Q = self.Wq_1(x)  # [batch_size, seq_len, hidden_dim]
        K = self.Wk_1(x)  # [batch_size, seq_len, hidden_dim]
        x_dp = torch.cat((F.dropout(gene_1, self.p_RNASeq), F.dropout(miRNA_1, self.p_miRNA)), 1)
        V = self.Wv_1(x_dp)  # [batch_size, seq_len, hidden_dim]

        # 求注意力分数
        attention_weights = torch.matmul(Q, Q.transpose(-2, -1)) / (50 ** 0.5)
        # 对注意力分数归一化，最后一个维度做
        attention_weights = nn.functional.softmax(attention_weights, dim=-1)
        # 将归一化后的权重与value相乘
        attended_values = torch.matmul(attention_weights, V)
        # print(attended_values)

        gene_2 = gene.unsqueeze(1)
        miRNA_2 = miRNA.unsqueeze(1)
        x = torch.cat((gene_2, miRNA_2), 1)

        Q_2 = self.Wq_2(x)  # [batch_size, seq_len, hidden_dim]
        K_2 = self.Wk_2(x)  # [batch_size, seq_len, hidden_dim]
        x_dp = torch.cat((F.dropout(gene_2, self.p_RNASeq), F.dropout(miRNA_2, self.p_miRNA)), 1)
        V_2 = self.Wv_2(x_dp)  # [batch_size, seq_len, hidden_dim]

        # 求注意力分数
        attention_weights_2 = torch.matmul(Q_2, Q_2.transpose(-2, -1)) / (50 ** 0.5)
        # 对注意力分数归一化，最后一个维度做
        attention_weights_2 = nn.functional.softmax(attention_weights_2, dim=-1)
        # 将归一化后的权重与value相乘
        attended_values_2 = torch.matmul(attention_weights_2, V_2)
        temp_1 = self.get_scale_feature(self.atten_fuse_weight(
            torch.cat((torch.squeeze(attended_values[:, 0, :]), torch.squeeze(attended_values[:, 1, :])), 1)))
        temp_2 = self.get_scale_feature(self.atten_fuse_weight(
            torch.cat((torch.squeeze(attended_values_2[:, 0, :]), torch.squeeze(attended_values_2[:, 1, :])), 1)))

        return self.BN(torch.squeeze(torch.cat((F.normalize(temp_1), F.normalize(temp_2)), 1)))

    def get_scale_feature(self, temp_fea):
        return torch.sqrt(F.relu(temp_fea)) - torch.sqrt(F.relu(-temp_fea))

    def get_survival_result(self, gene, miRNA, W, W_assist):
        gene_NN = self.get_gene_feature(gene)
        miRNA_NN = self.get_miRNA_feature(miRNA)
        gene_miRNA = self.get_attention_feature(gene_NN, miRNA_NN)
        final_fea = self.get_graph_feature(gene_miRNA, W, W_assist)
        return self.Cox_layer(final_fea)

    def get_attention_feature(self, RNASeq_feature, miRNA_feature):
        mul_head_feature = self.get_Transformer_feature(RNASeq_feature, miRNA_feature)
        return mul_head_feature

