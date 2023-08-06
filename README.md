# PMN_analysis
Codes for specific scRNA-seq analysis on PMN

'cell_crosstalk_enrichment.R' : Test if a receptor-ligand pair with significant signals frequently occurred in cell-cell communications either in PMN cells or controls using an enrichment method proposed by PMID: 34653364. The input files are nets/signaling pairs predicted to be significantly different between PMN and controls using CellChat(PMID:33597522).

'./PDM/calculate_PDM_for_podo.R' :Identify the potential tipping state of 114 podocytes from PMN patients. The PMN podocytes were subclustered into three states along PMN processing beforehand, and the control podocytes were used as the initial state (state 0) in the comparison. Gene modules were constructed by known regulons and were listed in MARAregulon.RData. To implement 'calculate_PDM_for_podo.R', several R packages, i.e. 'dplyr','tidyverse','MASS' and 'rrcov', should be installed.  
