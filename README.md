# CoFormerSurv
an innovative collaborative Transformer framework for multi-omics survival analysis

# Environment Installation

$ conda env create -f environment.yml

$ source activate survival_analysis

# Data Download

The TCGA datasets are publicly available and can be obtained from the UCSC Xena database (https://xena.ucsc.edu/). The original data can be accessed at: https://www.cancer.gov/ccg/research/genome-sequencing/tcga.


# Data Preprocessing

#### cd Federated-Transfer-MOSAHit/Data_preprocessing

The python program (extract_RNASeq_expression.py) summarized the gene expression of patients with the same type of cancer into a csv file. 

#### python extract_RNASeq_expression.py

#### python extract_miRNA_expression.py

The python program (RNASeq_variable_filter.py) imputed missing data, filtered out noise-sensitive features (those whose values ​​remained almost constant across samples), and performed a logarithmic transformation.

#### python RNASeq_variable_filter.py

#### python microRNA_variable_filter.py

The python program (clinical_preprocess.py)  aggregated the survival data of cancer patients (such as censoring indicator δ and observed time O).

#### python clinical_preprocess.py

The python program (matched_patients.py) only kept the patients with matched gene expression, microRNA expression data and survival data.

#### python matched_patients.py

# Run the main routine

#### git clone https://github.com/LiminLi-xjtu/Federated-Transfer-MOSAHit.git

#### cd Federated-Transfer-MOSAHit

#### python DPFL.py
