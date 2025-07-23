# -*- coding: utf-8 -*-
"""
Created on Wed Jul 23 12:55:28 2025

The goal of this script to take in data from Zhuang et al. (see other files) 
and be able to explore spatial transcriptomics of AOB MCs, specifically. 

This will be done in a few steps including:
 --> data loading
 --> joining dataframes
 --> filtering for excitatory neurons within the accessory olfactory bulb
 --> visualizing genetic diversity of these neurons

@author: Kevin
"""
#%% imports
import pandas as pd
from pathlib import Path
import numpy as np
import anndata
import time
import matplotlib.pyplot as plt

from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache
from AllenBrainAtlasPlottingFunctions import plot_sections, plot_heatmap
from combining_filtering_abc_data import create_cell_extended, filter_cell_extended, add_gene_expression

#%% setting up info for data
# labels for all Zhuang datasets
datasets = ['Zhuang-ABCA-1', 'Zhuang-ABCA-2', 'Zhuang-ABCA-3']
# loading / caching relevant data:
download_base = Path('C:/Users/IGD/Documents/KJM/data/abc_10x/data/abc_atlas')
abc_cache = AbcProjectCache.from_cache_dir(download_base)

#%% combine datasets with other metadata
cell_extended = create_cell_extended(datasets)

#%% filter for only AOB glutamatergic cells
cell_extended_filt = filter_cell_extended(cell_extended,'AOB','Glut');

#%% combining cell_extended dataframes with gene expression
# loading in gene data
genes = abc_cache.get_metadata_dataframe(directory=datasets[0],
                                        file_name='gene')
genes.set_index('gene_identifier', inplace=True)

# you can filter the genes based on a set list, but it must be a dataframe
#for example, here we are selecting a subset of screened genes
gnames = ['Slc17a6','Slc32a1','Gnrh1','Esr1','Esr2','Npy1','Npy2','Cyp19a1','Tac2','Crhr1','Crhr2']
gene_bool = [x in gnames for x in genes.gene_symbol]
genes_filtered = genes[gene_bool]

cell_extended_filt_genex = add_gene_expression(cell_extended_filt,genes_filtered)

#%% with all data compiled, we can combine across dataframes as follows
df_all_datasets = pd.DataFrame()
for k in cell_extended_filt_genex:
    df_all_datasets = pd.concat([df_all_datasets,
                                 cell_extended_filt_genex[k]])
    
    
#%% data visualization

        