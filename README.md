# PAL description

Pathway Analysis of Longitudinal data (PAL) is a pathway analysis tool that provides pathway values for each analysed sample and pathway separately. PAL allows the analysis of complex study designs such as longitudinal data. Details about the algorithm are available in the original open access publication "Pathway analysis of longitudinal data with case studies in early type 1 diabetes" [1] (please cite it if you utilise PAL in your research). In case of bugs, missing documentation, or problems with installation, please contact us: maria.jaakkola@utu.fi

PAL can be installed by opening R and typing devtools::install_github("elolab/PAL") (requires package devtools to be installed). Notably, usage of PAL requires installation of R package PASI (devtools::install_github("elolab/PASI"))

### Input and output

The only mandatory input from the user are **data** and **grouplabels**.
| Input | Description |
| ----------- | ----------- |
| data | A data frame of normalized gene expression data in log2 scale (rows: Entrez genes, cols: samples) |
| grouplabels | An integer vector indicating sample groups. Control samples should be indicated with label 0. |
| pathwayadress | A directory path to the user's own pathways to be analysed on top of the KEGG pathways. If NULL (default), only the automatically accessed KEGG pathways are analysed. |
| datatype | Either "microarray" or "rnaseq" (default). This is used only if argument **noisedefault** is set to "automatic". |
| noisedefault | The default cutoff (numeric value or character ”automatic” (default)) for real signal. If set to NA, no lowly expressed genes are filtered out. |
| score | Either "activity" (default) or "deregulation" based on what the returned pathway scores should reflect (more details below). |
| nodemin | Indicates the minimum nuber of measured nodes in a pathway to be analysed (default 5). Pathways with fewer measured nodes are excluded from the analysis. |
| info | A data frame including all variables to be used in the model fitting in the neutralisation step as columns. Rows correspond to samples in **data**. Rows and columns should be named. |
| neutralize | A logical vector indicating which coefficients in 'info' should be neutralised (T=neutralise, F=use in model fitting, but don't neutralise). The length should match the number of columns in **info**. |
| mainfeature | A numeric or character vector corresponding to samples (cols in 'data'). If provided, pathways' significance levels according to this coefficient are returned.|

In case the expression data contains no sample groups, **grouplabels** can be set to dummy value of only zeros rep(0,ncol(data)).

There are two options for argument **score**. The default one is "activity" and if it is selected, the final pathway scores indicate how active each pathway is in comparison to the other samples. Negative value indicate inactivity, value close to 0 normal activity, and positive value high activity. If **score** is set to "deregulation", the pathway scores indicate how normally the pathways behave as compared to a typical control sample. Value close to 0 means that the pathway is not deregulated and a high value means that it is.

The main function PAL returns a list with two or three elements: a matrix of pathway scores, a matrix of significance levels and coefficients (only is argument **mainfeature** is provided) of the analysed pathways, and an info matrix describing the analysed pathways. The pathway scores can be used for further analyses.


##### References

[1] Add publication info when known


