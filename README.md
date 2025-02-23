# Overview of Incytr
Signaling pathways discovery and quantitative analysis from multi-modal data
(Available on biorxiv: https://www.biorxiv.org/content/10.1101/2025.02.06.636961v1)

Incytr is an open-source R package to infer cell signaling pathways between cell types with the structure "Ligand-Receptor-Effector Molecule-Target" using multi-modal data and prior knowledge. 

Incytr identifies cell signaling pathways from scRNA-seq data alone or integrated with proteomics, phosphoproteomics, and kinase-substrate specificity. The required inputs for Incytr are single-cell transcriptomics data, cell group labels, and user-selected sender and receiver genes which can be any genes measured in the data. The optional inputs include condition labels for the cells, proteomics, phosphoproteomics data, and a predicted kinase-substrate list (see Kinase-Substrate Matching in Supplementary Materials). Below is an overview of the Incytr workflow:

<img width="536" alt="image" src="https://github.com/user-attachments/assets/2c20a90c-7c73-4ef4-b4ea-346ce70910b9" />


