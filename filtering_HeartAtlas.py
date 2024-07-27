import scanpy as sc

scanpy_obj=sc.read_h5ad("~/MSc/project/deconvolution/global_files/Global_raw.h5ad")
scanpy_obj.obs
if 'Celltype' in scanpy_obj.obs:
    print('"Celltype" is present in the obs attribute.')
else:
    print('"Celltype" is not present in the obs attribute.')
    
#QC checks
#remove cells with less than 200 genes and genes expressed in less than 3 cells
sc.pp.filter_cells(scanpy_obj,min_genes=200)
sc.pp.filter_genes(scanpy_obj,min_cells=3)
#normalising the data
sc.pp.normalize_total(scanpy_obj,target_sum=1e4)
sc.pp.log1p(scanpy_obj)

#renaming to fit obs specifications of Scaden
if 'Celltype' not in scanpy_obj.obs.columns:
    raise ValueError('"Celltype" is not present in the obs attribute.')

print(scanpy_obj.obs.head())

old_name = 'cell_type'
new_name = 'Celltype'
if old_name in scanpy_obj.obs.columns:
    scanpy_obj.obs.rename(columns={old_name: new_name}, inplace=True)
else:
    print(f"Column '{old_name}' not found in obs.")

#identifying code for left ventricle
print(scanpy_obj.obs['region'].unique())

import numpy as np
from scipy.sparse import issparse

def sample_scaling(scanpy_obj, scaling_option):
    if issparse(scanpy_obj):
        scanpy_obj = scanpy_obj.toarray()
        scanpy_obj = np.log2(scanpy_obj + 1)
    else:
        scanpy_obj = np.log2(x + 1)
    return scanpy_obj
region_to_filter = 'LV'
filtered_data = scanpy_obj[scanpy_obj.obs['region'] == region_to_filter]
filtered_data.write_h5ad('LV_py.h5ad')
