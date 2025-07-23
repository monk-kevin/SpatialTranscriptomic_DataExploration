 """ 
This set of functions takes publically-available datasets and processes them 
for additional analysis. Specifically:
    
    1) <create_cell_extended> combines across different files to incorporate 
    cellular identity along with relevant metadata
    
    2) <filter_cell_extended> filters larger dataframes to select specific cell
    populations defined by the neurotransmitter used and location in the brain.
    
    3) <add_gene_expression> combines a dataframe with cellular identity with 
    known transcriptomic expression pulled in from a different data stream.
    
All data is provided by the Allen Brain Institute. Code written here has been
based off of code they provide to suit our custom needs.
"""

import pandas as pd
from pathlib import Path
import numpy as np
import anndata
import time
import matplotlib.pyplot as plt

from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

def create_cell_extended(datasets):
    """
    This function loads and creates an extended dataframe that compiles cell identity
    with a variety of other factors including:
        cluster id
        cluster color
        ABC parcellated structure
    """
    # loading / caching relevant data:
    download_base = Path('C:/Users/IGD/Documents/KJM/data/abc_10x/data/abc_atlas')
    abc_cache = AbcProjectCache.from_cache_dir(download_base)
    print(f'Current Manifest is: {abc_cache.current_manifest}')
        
    #reading in cluster details and providing color for consistent labeling
    cluster_details = abc_cache.get_metadata_dataframe(
        directory='WMB-taxonomy',
        file_name='cluster_to_cluster_annotation_membership_pivoted',
        keep_default_na=False
        )
    cluster_details.set_index('cluster_alias', inplace=True)

    cluster_colors = abc_cache.get_metadata_dataframe(
        directory='WMB-taxonomy',
        file_name='cluster_to_cluster_annotation_membership_color',
        )
    cluster_colors.set_index('cluster_alias', inplace=True)


# loading in cell metadata for each dataframe (a dataframe given by the datasets above)
    print('Creating components...')
    cell = {}
    ccf_coordinates = {}
    for d in datasets :
        #for each dataset, 
        # load in the metadata
        cell[d] = abc_cache.get_metadata_dataframe(
            directory=d,
            file_name='cell_metadata',
            dtype={"cell_label": str}
        )
        cell[d].set_index('cell_label', inplace=True)
            
        sdf = cell[d].groupby('brain_section_label')
            
        print(d,":","Number of cells = ", len(cell[d]), ", ", "Number of sections =", len(sdf))
            
        # load in the coordinate information
        ccf_coordinates[d] = abc_cache.get_metadata_dataframe(directory=f"{d}-CCF", file_name='ccf_coordinates')
        ccf_coordinates[d].set_index('cell_label', inplace=True)
        ccf_coordinates[d].rename(columns={'x': 'x_ccf',
                                           'y': 'y_ccf',
                                           'z': 'z_ccf'},
                                  inplace=True)
        
    # parcellation data
    parcellation_annotation = abc_cache.get_metadata_dataframe(directory="Allen-CCF-2020",
                                                               file_name='parcellation_to_parcellation_term_membership_acronym')
    parcellation_annotation.set_index('parcellation_index', inplace=True)
    parcellation_annotation.columns = ['parcellation_%s'% x for x in  parcellation_annotation.columns]

    parcellation_color = abc_cache.get_metadata_dataframe(directory="Allen-CCF-2020",
                                                              file_name='parcellation_to_parcellation_term_membership_color')
    parcellation_color.set_index('parcellation_index', inplace=True)
    parcellation_color.columns = ['parcellation_%s'% x for x in  parcellation_color.columns]
        
    
    print('Creating cell_extended...')
    # joining previous dataframes with cluster details and colors
    cell_extended = {}
    for d in datasets :
        cell_extended[d] = cell[d].join(cluster_details, on='cluster_alias')
        cell_extended[d] = cell_extended[d].join(cluster_colors, on='cluster_alias')
        cell_extended[d] = cell_extended[d].join(ccf_coordinates[d], how='inner')
        cell_extended[d] = cell_extended[d].join(parcellation_annotation, on='parcellation_index')
        cell_extended[d] = cell_extended[d].join(parcellation_color, on='parcellation_index')
            
    return cell_extended

def filter_cell_extended(cell_extended,
                         structure_mask = 'AOB',
                         nt_mask = 'Glut'):
    """
    This function filters the large extended data frame to only consider 
    neurons within the defined area (set by structure_mask, AOB is default) and
    by neurotransmitter used (set by nt_mask, Glut by default)
    """
    
    # creating a new dataframe only with cells within the olfactory bulb
    cell_extended_filt = {}

    for key in cell_extended:
        curr_df = cell_extended[key]
        # struct_mask_bin = curr_df['parcellation_structure'] == structure_mask
        # nt_mask_bin = curr_df['neurotransmitter'] == nt_mask
        # combined_mask = struct_mask_bin + nt_mask_bin;
        
        # cell_extended_filt[key] = curr_df.iloc[list(combined_mask)]
        
        cell_extended_filt[key] = curr_df.loc[
            (curr_df['parcellation_structure'] == structure_mask) &
            (curr_df['neurotransmitter'] == nt_mask)]

        print(f'{key}: {cell_extended_filt[key].shape[0]} {nt_mask}-ergic neurons within {structure_mask}')
        
    return cell_extended_filt
    
def add_gene_expression(cell_extended,genes):
    """
    This function incorporates genetic expression data into 
    a cell_extended dictionary of dataframes. cell_extended is a dictionary of 
    dataframes and genes is a dataframe with information for screened genes.
    """
    
    print(f'Finding expression data for {len(genes)} genes')
    download_base = Path('C:/Users/IGD/Documents/KJM/data/abc_10x/data/abc_atlas')
    abc_cache = AbcProjectCache.from_cache_dir(download_base)
    
    cell_w_genex = {}
    for key in cell_extended:    
        file = abc_cache.get_data_path(directory=key, 
                                       file_name=f"{key}/log2")
        
        adata = anndata.read_h5ad(file, backed='r')
        
        start = time.process_time()
        gdata = adata[:, genes.index].to_df()
        gdata.columns = genes.gene_symbol
        cell_w_genex[key] = cell_extended[key].join(gdata)
        
        print(key,"-","time taken: ", time.process_time() - start)
        
        adata.file.close()
        del adata
    
    return cell_w_genex