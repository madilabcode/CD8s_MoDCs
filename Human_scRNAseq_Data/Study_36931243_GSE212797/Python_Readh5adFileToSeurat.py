# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 11:44:38 2024

@author: kfiri
"""
#%%
import os
path = os.path.dirname(os.path.abspath(__file__))
os.chdir(path)
#%%
import scanpy as sc
from scipy import io
import gzip

#file = '.\GSE212797_adata.h5ad'
file = '.\GSE212797_adata_CD8.h5ad'
adata = sc.read_h5ad(file)
adata = adata.raw.to_adata()  #only if adata has RAW saved and thats what you want!!
adata.obs_names
adata.var_names
adata.obs
dir_to_save = r'.\GSE212797_adata_CD8'
#%%
#barcodes
with open(dir_to_save + r'\barcodes.tsv', 'w') as f:
    for item in adata.obs_names:
        f.write(item + '\n')
#%%
#features
with open(dir_to_save + r'\features.tsv', 'w') as f:
    for item in ['\t'.join([x,x,'Gene Expression']) for x in adata.var_names]:
        f.write(item + '\n')
#%%
#matrix
io.mmwrite(dir_to_save + r'\matrix', adata.X.T)
#%%

# Get a list of all files in the directory
files = os.listdir(dir_to_save)

# Loop through each file and compress it using gzip
for file_name in files:
    file_path = os.path.join(dir_to_save, file_name)
    with open(file_path, 'rb') as f_in:
        with gzip.open(file_path + '.gz', 'wb') as f_out:
            f_out.writelines(f_in)
#%%
#annotations
adata.obs.to_csv(dir_to_save + r'\metadata.csv')

