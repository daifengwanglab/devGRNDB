# Converting h5ad pseudocells files to loom files for pyscenic runs.. The h5ad files contain umi counts.

import numpy as np
import scanpy as sc
import pandas as pd
import loompy as lp

# case information
type = 'Telencephalon'
linType = 'LGE'
region = 'subtypes'
geneList = pd.read_csv('GeneList.csv') # the consensus gene list used for GRN construction


initPath = 'data/subregions/'+type+'_Pcells.h5ad'
metaAvail = True
savePath = 'human/'+region+'/'+type+'/'+type+'_'+linType
metaPath = 'data/meta/Harmony_pr_pseudotime_rep4_'+type+'_'+linType+'_v3_01272025.csv'
adata = sc.read_h5ad(initPath)
adata.obs_names = adata.obs['cell_ID']
print(adata.obs_names[0:10])
print(adata.shape)



sc.pp.normalize_total(adata, target_sum = 1e6)

print(adata.shape)
if metaAvail == True:
	meta = pd.read_csv(metaPath,index_col = 'Unnamed: 0')
	print(meta.shape)
	adata = adata[adata.obs_names.isin(meta.index)]

adata.write_h5ad(savePath+'.h5ad')

adata = adata[:, adata.var_names.isin(geneList['x'])]

row_attrs = { 
    "Gene": np.array(adata.var.index) ,
}
col_attrs = { 
    "CellID":  np.array(adata.obs.index) ,
    "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
}

lp.create(savePath+'.loom', adata.X.transpose(), row_attrs, col_attrs )
