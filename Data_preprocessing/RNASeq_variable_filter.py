import os
import pandas as pd
import numpy as np
# 导入CSV安装包
# import csv



cancer_list = ['HNSC']
source_file_path = 'Collaborative_Transformer/Datasets_and_Preprocessing/data/RNASeq_source_to_csv'
target_file_path = 'Collaborative_Transformer/Datasets_and_Preprocessing/data/RNASeq_variable_filter_csv/'


for root, dirs, files in os.walk(source_file_path):
    for file in files:
        cancer_type = file.split('.')[0].split('-')[1]

        if cancer_type in cancer_list:
            file_path = str(os.path.join(root, file).encode('utf-8'), 'utf-8')
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
                if col[:2] == 'EN':
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
                if col[:2] == 'EN':
                    with open('gene_process' + '.log', 'a') as f:
                        f.writelines('gene_name:' + col + '\n')
                    temp[col] = temp[col].astype(float)
                    log2_val = np.log2(np.array(temp[col]).squeeze() +1)
                    var = np.var(log2_val)
                    temp[col] = log2_val.tolist()
                    var_value_list.append(var)
                    colname_list.append(col)
            ddff = pd.DataFrame({'col_name':colname_list, 'var_value': var_value_list})
            var_col_name = np.squeeze(ddff.sort_values(by='var_value', ascending=False).iloc[:10000][['col_name']].values).tolist()
            print(var_col_name)
            temp['submitter_id'] = temp['submitter_id'].str[0:12]
            temp = temp.drop_duplicates('submitter_id')
            var_col_name.append('label')
            # var_col_name.append('file_id')
            # var_col_name.append('file_name')
            var_col_name.append('submitter_id')
            print(temp['label'])
            print(cancer_type + "筛选后gene的数目：" + str(len(var_col_name)-2))
            RNASeq_expresssion = temp[var_col_name]
            RNASeq_expresssion[RNASeq_expresssion['label'].isin([file.split('.')[0]])].to_csv(target_file_path + file, index=False, header=True)
