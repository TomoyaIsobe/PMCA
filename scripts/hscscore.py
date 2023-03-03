import pandas as pd
import numpy as np
import pickle
import sklearn
import platform
import scanpy as sc
import warnings
warnings.filterwarnings("ignore")

# params
input_adata = 'data/PMCA_Jak2_raw_count.h5ad'
outdir = 'data/'
outname_prefix = 'PMCA_Jak2'
pkl_file = 'scripts/hscScore_model.pkl'
model_genes_file = 'scripts/model_molo_genes.txt'

# loading data
adata = sc.read(input_adata)

# loading model files
hsc_score = pickle.load(open(pkl_file, 'rb'))
model_genes = np.genfromtxt(model_genes_file, dtype='str')

# preprocessing
OLgenes = np.intersect1d(model_genes, adata.var_names)
adata_sub = adata[:,OLgenes].copy()
count_data_molo = pd.DataFrame(adata_sub.X.todense(), index=adata_sub.obs_names, columns=adata_sub.var_names)

missingGenes = np.setdiff1d(model_genes, adata.var_names)
for g in missingGenes:
    count_data_molo[g] = [0]*count_data_molo.shape[0]

count_data_molo = count_data_molo.loc[:,model_genes].copy()

def total_count_normalise(count_matrix):
    """Normalise count matrix for input into hscScore model.
    Performs read depth normalisation normalising each cell so that normalised 
    counts sum to the same value.
    
    Parameters
    ----------
    count_matrix : pandas dataframe
        Gene count matrix of dimension cells x genes with column names as genes
        and index as cell names
    
    Returns
    -------
    **norm_matrix** : pandas dataframe
        Normalised count matrix of dimension cells x genes
    """
    
    # Set the value normalised counts will sum to for each cell
    wilson_molo_genes_median_counts = 18704.5
    
    # Scale rows
    count_matrix_expression = np.array(count_matrix, dtype='float')
    counts_per_cell = np.sum(count_matrix_expression, axis=1)
    counts_per_cell += (counts_per_cell == 0)
    counts_per_cell /= wilson_molo_genes_median_counts
    norm_matrix_expression =  count_matrix_expression/counts_per_cell[:, None]
    norm_matrix = pd.DataFrame(norm_matrix_expression, index=count_matrix.index,
                               columns=count_matrix.columns)
    # log + 1 transform the data
    norm_matrix = np.log(norm_matrix + 1)
    
    return norm_matrix


normalised_data_molo = total_count_normalise(count_data_molo)
predicted_hsc_scores = hsc_score.predict(np.array(normalised_data_molo))

predicted_hsc_scores_df = pd.DataFrame(predicted_hsc_scores)
predicted_hsc_scores_df.index = list(adata.obs.index)
predicted_hsc_scores_df.to_csv(outdir+outname_prefix+'_hsc_scores.csv')

print('Completed')
