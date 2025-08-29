import os
import pandas as pd
import numpy as np
# 导入CSV安装包
# import csv

# ['KIRC', 'KIRP','KICH']
# cancer_type = 'Kindey'
# base_path = '/share/home/4120107034/Collaborative_Transformer/Datasets_and_Preprocessing/data'
# file_path1 = base_path + '/miRNA_source_to_csv/TCGA-KIRC.csv'
# file_path2 = base_path + '/miRNA_source_to_csv/TCGA-KIRP.csv'
# file_path3 = base_path + '/miRNA_source_to_csv/TCGA-KICH.csv'
# target_file_path = base_path + '/miRNA_variable_filter_csv'
#
# temp1 = pd.read_csv(file_path1, sep=',')
# temp2 = pd.read_csv(file_path2, sep=',')
# temp3 = pd.read_csv(file_path2, sep=',')
# temp = pd.concat([temp1, temp2, temp3], axis=0)
#
# print(cancer_type + " gene_feature的数目：" + str(temp.shape[1]))
# # temp = temp.dropna(axis=0, how='any')
# temp = temp.dropna(axis=1, thresh=temp.shape[0] * 0.9)
# print("删除None后" + cancer_type + " gene_feature的数目：" + str(temp.shape[1]))
#
# print(cancer_type + "样本的数目：" + str(temp.shape[0]))
# # temp = temp.dropna(axis=0, how='any')
# temp = temp.dropna(axis=0, thresh=temp.shape[1] * 0.5)
# print("删除None后" + cancer_type + "样本的数目：" + str(temp.shape[0]))
# columns_name = list(temp.columns)
# for col in columns_name:
#     if col[:3] == 'hsa':
#         temp[col] = np.squeeze(np.array(temp[[col]].fillna(temp[col].median())))
# # if flag:
# #     temp = temp_expression
# #     flag = False
# # else:
# #     temp = temp.append(temp_expression)
# # break
#
# colname_list = []
# var_value_list = []
# for col in columns_name:
#     if col[:3] == 'hsa':
#         temp[col] = temp[col].astype(float)
#         var = np.var(np.log2(np.array(temp[col]).squeeze() + 1))
#         temp[col] = np.log2(np.array(temp[col]).squeeze() + 1).tolist()
#         var_value_list.append(var)
#         colname_list.append(col)
# ddff = pd.DataFrame({'col_name':colname_list, 'var_value': var_value_list})
# var_col_name = np.squeeze(ddff.sort_values(by='var_value', ascending=False).iloc[:600][['col_name']].values).tolist()
# print(var_col_name)
# temp['submitter_id'] = temp['submitter_id'].str[0:12]
# temp = temp.drop_duplicates('submitter_id')
# var_col_name.append('label')
# # var_col_name.append('file_id')
# # var_col_name.append('file_name')
# var_col_name.append('submitter_id')
# print("筛选后miRNA的数目：" + str(len(var_col_name)-2))
# miRNA_expresssion = temp[var_col_name]
# miRNA_expresssion.to_csv(target_file_path + '/Kindey.csv', index=False, header=True)
#
# # ['LUSC', 'LUAD']
# cancer_type = 'Lung'
# base_path = '/share/home/4120107034/Collaborative_Transformer/Datasets_and_Preprocessing/data'
# file_path1 = base_path + '/miRNA_source_to_csv/TCGA-LUSC.csv'
# file_path2 =base_path + '/miRNA_source_to_csv/TCGA-LUAD.csv'
# target_file_path = base_path + '/miRNA_variable_filter_csv'
#
# temp1 = pd.read_csv(file_path1, sep=',')
# temp2 = pd.read_csv(file_path2, sep=',')
# temp = pd.concat([temp1, temp2], axis=0)
#
# print(cancer_type + " gene_feature的数目：" + str(temp.shape[1]))
# # temp = temp.dropna(axis=0, how='any')
# temp = temp.dropna(axis=1, thresh=temp.shape[0] * 0.9)
# print("删除None后" + cancer_type + " gene_feature的数目：" + str(temp.shape[1]))
#
# print(cancer_type + "样本的数目：" + str(temp.shape[0]))
# # temp = temp.dropna(axis=0, how='any')
# temp = temp.dropna(axis=0, thresh=temp.shape[1] * 0.5)
# print("删除None后" + cancer_type + "样本的数目：" + str(temp.shape[0]))
# columns_name = list(temp.columns)
# for col in columns_name:
#     if col[:3] == 'hsa':
#         temp[col] = np.squeeze(np.array(temp[[col]].fillna(temp[col].median())))
# # if flag:
# #     temp = temp_expression
# #     flag = False
# # else:
# #     temp = temp.append(temp_expression)
# # break
#
# colname_list = []
# var_value_list = []
# for col in columns_name:
#     if col[:3] == 'hsa':
#         temp[col] = temp[col].astype(float)
#         var = np.var(np.log2(np.array(temp[col]).squeeze() + 1))
#         temp[col] = np.log2(np.array(temp[col]).squeeze() + 1).tolist()
#         var_value_list.append(var)
#         colname_list.append(col)
# ddff = pd.DataFrame({'col_name':colname_list, 'var_value': var_value_list})
# var_col_name = np.squeeze(ddff.sort_values(by='var_value', ascending=False).iloc[:600][['col_name']].values).tolist()
# print(var_col_name)
# temp['submitter_id'] = temp['submitter_id'].str[0:12]
# temp = temp.drop_duplicates('submitter_id')
# var_col_name.append('label')
# # var_col_name.append('file_id')
# # var_col_name.append('file_name')
# var_col_name.append('submitter_id')
# print("筛选后miRNA的数目：" + str(len(var_col_name)-2))
# miRNA_expresssion = temp[var_col_name]
# miRNA_expresssion.to_csv(target_file_path + '/Lung.csv', index=False, header=True)


# flag = True
# cancer_list = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KIRC', 'KIRP', 'KICH',
#                'LAML', 'LGG', 'LIHC', 'LUSC', 'LUAD', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM',
#                'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']
cancer_list = ['HNSC']
source_file_path = '/share/home/4120107034/Collaborative_Transformer/Datasets_and_Preprocessing/data/miRNA_source_to_csv'
target_file_path = '/share/home/4120107034/Collaborative_Transformer/Datasets_and_Preprocessing/data/miRNA_variable_filter_csv/'
for root, dirs, files in os.walk(source_file_path):
    for file in files:
        cancer_type = file.split('.')[0].split('-')[1]
        if cancer_type in cancer_list:
            file_path = str(os.path.join(root, file).encode('utf-8'), 'utf-8')
            print(file_path)
            temp = pd.read_csv(file_path, sep=',')
            print(cancer_type + " gene_feature的数目：" + str(temp.shape[1]))
            # temp = temp.dropna(axis=0, how='any')
            temp = temp.dropna(axis=1, thresh=temp.shape[0] * 0.9)
            print("删除None后" + cancer_type + " gene_feature的数目：" + str(temp.shape[1]))

            print(cancer_type + "样本的数目：" + str(temp.shape[0]))
            # temp = temp.dropna(axis=0, how='any')
            temp = temp.dropna(axis=0, thresh=temp.shape[1] * 0.5)
            print("删除None后" + cancer_type + "样本的数目：" + str(temp.shape[0]))
            columns_name = list(temp.columns)
            for col in columns_name:
                if col[:3] == 'hsa':
                    temp[col] = np.squeeze(np.array(temp[[col]].fillna(temp[col].median())))
            # if flag:
            #     temp = temp_expression
            #     flag = False
            # else:
            #     temp = temp.append(temp_expression)
            # break

            colname_list = []
            var_value_list = []
            for col in columns_name:
                if col[:3] == 'hsa':
                    temp[col] = temp[col].astype(float)
                    var = np.var(np.log2(np.array(temp[col]).squeeze() + 1))
                    temp[col] = np.log2(np.array(temp[col]).squeeze() + 1).tolist()
                    var_value_list.append(var)
                    colname_list.append(col)
            ddff = pd.DataFrame({'col_name':colname_list, 'var_value': var_value_list})
            var_col_name = np.squeeze(ddff.sort_values(by='var_value', ascending=False).iloc[:600][['col_name']].values).tolist()
            print(var_col_name)
            temp['submitter_id'] = temp['submitter_id'].str[0:12]
            temp = temp.drop_duplicates('submitter_id')
            var_col_name.append('label')
            # var_col_name.append('file_id')
            # var_col_name.append('file_name')
            var_col_name.append('submitter_id')
            print("筛选后miRNA的数目：" + str(len(var_col_name)-2))
            miRNA_expresssion = temp[var_col_name]
            miRNA_expresssion[miRNA_expresssion['label'].isin([file.split('.')[0]])].to_csv(target_file_path + file,
                                                                                     index=False, header=True)