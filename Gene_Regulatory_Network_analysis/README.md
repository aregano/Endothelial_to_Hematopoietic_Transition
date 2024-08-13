# Gene Regulatory Network Analysis with CellOracle

In order to check for alterations in the gene networks of the AGM upon the various conditions I conducted GRN analysis using CellOracle (https://morris-lab.github.io/CellOracle.documentation/)

Briefly I selected the cells that belonged to the AGM and extracted GRN links using a publicly available scATAC dataset that was performed on the AGM (GSE137115 accession number from Zhu et al. (2020))

I looked into the differences in GRN between conditions by performing Network Analysis while looking at the cell populations withing the AGM. (Network_analysis_scAGM_EHT.ipynb)

Furthermore I performed in silico TF perturbation analysis, where CellOracle uses the GRN models to simulate cell identity shifts in response to TF perturbation. I performed the Perturbation analysis with various genes selected for a count of 0: _Hey1, Hey2, Hes1, Gata2, Runx1, Mycn, Myc_ (KO_Perturbations.Hey1_KO_simulation_with_scAGM_EHT.ipynb)
