import os
import os.path
import gzip
import numpy as np
import pandas as pd


def read_gz_file(path):
    if os.path.exists(path):
        with gzip.open(path, 'r') as pf:
            for line in pf:
                yield line
    else:
        print('the path [{}] is not exist!'.format(path))


# # 导入CSV安装包
# import csv
# 导入json包
import json

# # 创建文件对象
# file_write_name = open('tcga_table_acc.csv', 'w', newline='', encoding='utf-8')
# # 基于文件对象构建 csv写入对象
# csv_writer = csv.writer(file_write_name)

path = []
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-ACC")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-BLCA")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-BRCA")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-CESC")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-CHOL")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-COAD")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-DLBC")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-ESCA")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-GBM")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-HNSC")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-KICH")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-KIRC")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-KIRP")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-LAML")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-LGG")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-LIHC")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-LUAD")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-LUSC")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-MESO")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-OV")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-PAAD")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-PCPG")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-PRAD")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-READ")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-SARC")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-SKCM")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-STAD")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-TGCT")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-THCA")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-THYM")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-UCEC")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-UCS")
path.append("/home/datasets/TCGA/RNA-Seq/TCGA-UVM")

target_path = './Datasets_and_Preprocessing/data/RNASeq_source_to_csv/'

for type_path in path:
    label = type_path.split('/')[5]
    print(label)
    # metadata = open(path + "/metadata.cart.2021-10-09.json", encoding='utf-8')
    metadata = open(type_path + "/metadata.cart.2021-10-09.json", encoding='utf-8')
    json_data = json.load(metadata)
    # print(json_data)
    file_name_list = []
    file_id_list = []
    submitter_id_list = []
    for list_data in json_data:
        # print(list_data['file_name'])
        file_name = list_data['file_name']
        if file_name is None:
            file_name = "NA"
        file_name_list.append(file_name)
        file_id = list_data['file_id']
        if file_id is None:
            file_id = "NA"
        file_id_list.append(file_id)
        submitter_id = list_data['associated_entities'][0]['entity_submitter_id']
        if submitter_id is None:
            submitter_id = "NA"
        submitter_id_list.append(submitter_id)
    data = {"file_name": file_name_list, "file_id": file_id_list, "submitter_id": submitter_id_list}
    meta_dataframe = pd.DataFrame(data)
    # for root, dirs, files in os.walk(path + "/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification"):
    count = 0
    flag = True
    for root, dirs, files in os.walk(type_path + "/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification"):
        for file in files:
            # print(str(os.path.join(root, file).encode('utf-8'), 'utf-8')[-2:])
            file_path = str(os.path.join(root, file).encode('utf-8'), 'utf-8')
            col_list_name = []
            col_list_value = []
            if file_path[-2:] == 'gz':
                # print(file_path)
                file_id = file_path.split('/')[9]
                # print(file_id)
                file_name = file_path.split('/')[10]
                # print(file_name)
                if flag:
                    con = read_gz_file(file_path)
                    for line in con:
                        col_name = str(line, 'utf-8').split('\t')[0].split('.')[0]
                        # print(col_name)
                        if col_name[:2] == 'EN':
                            col_list_name.append(col_name)
                    col_list_name.append('label')
                    col_list_name.append('file_id')
                    col_list_name.append('file_name')
                    # csv_writer.writerow(col_list_name)
                    gene_dataframe = pd.DataFrame(columns=col_list_name)
                    flag = False
                con = read_gz_file(file_path)
                if getattr(con, '__iter__', None):
                    for line in con:
                        col_name = str(line, 'utf-8').split('\t')[0].split('.')[0]
                        if col_name[:2] == 'EN':
                            col_value = str(line, 'utf-8').split('\t')[1].replace('\n', '')
                            if col_value is None:
                                col_value = "NA"
                            # print(col_value)
                            col_list_value.append(col_value)
                    col_list_value.append(label)
                    col_list_value.append(file_id)
                    col_list_value.append(file_name)
                    # csv_writer.writerow(col_list_value)
                    count = count + 1
                    gene_dataframe.loc[count] = col_list_value
    print(label + ":" + str(count))
    #             print(gene_dataframe)
    # if count == 5:
    #     break
    # print(count)
    result = pd.merge(gene_dataframe, meta_dataframe, how="inner", on=['file_id', 'file_name'])
    result.to_csv(target_path + label + ".csv")

