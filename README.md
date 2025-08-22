# SpatialTranscriptomic_DataExploration
Code to investigate allen brain spatial transcriptomic data provided by the Allen Brain Institute (https://alleninstitute.github.io/abc_atlas_access/descriptions/Zhuang_dataset.html).

Allen Institue code can be found here (https://alleninstitute.github.io/abc_atlas_access/notebooks/zhuang_merfish_tutorial.html), some of which has been copied into the AllenBrainAtlasPlottingFunctions script within this repository.

AOB_MC_data_exploration is the main script to perform data processing and a first-pass of data visualization. To accomplish, I have created custom functions and saved them within the combining_filtering_abc_data.py file.

Using a combination of principal component analysis (PCA) and the k-means clustering algorithm, I am able to define at least 4 subtypes of AOB excitatory neurons:

<img width="826" height="819" alt="image" src="https://github.com/user-attachments/assets/19273d3b-ae8a-4a28-bfcd-6c946cf7c0ee" />

With these algorithms, I then performed feature importance to identify genes that are relevant for these genetic subtypes. Of particular note are Sox10, Reelin, and Cadherin-1.

Future work will be to determine the optimal set of clusters (e.g., varying the defined value, k) and visualize the reduced dataset by UMAP-embedding.
