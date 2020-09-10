# DEG
The experiment screened the differentially expressed genes of ALS and did some visual analysis.
Experimental materials include: human motor neuron tissue, mouse motor neuron, mouse oculomotor nucleus, mouse spinal cord tissue,which is used in GSE60856,GSE3343,GSE18579,GSE20598.

Take files in "homo" as an example:
GSE20589_series_matrix.txt is the profile data downloaded from GEO.
files in "result" includes DEGs filtered by FC and p-value,boxplot and volcano plot on DEGs,and pathways enriched from KEGG using DAVID.

To father use disease's DEGs and KS to calculate CMap score of compounds in CMap or LINCS ,see https://github.com/FlyInsect/pathway-and-lincs-data-KS-test
