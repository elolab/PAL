# PAL description

Pathway Analysis of Longitudinal data (PAL) is a pathway analysis tool that provides pathway values for each analysed sample and pathway separately. PAL allows the analysis of complex study designs such as longitudinal data and the effect of some variables (e.g. age) can be neutralised (see arguments **info** and **neutralize**) prior to the pathway analysis. Details about the algorithm are available in the original open access publication "Pathway analysis of longitudinal data with case studies in early type 1 diabetes" [1] (please cite it if you utilise PAL in your research). In case of bugs, missing documentation, or problems with installation, please contact us: maria.jaakkola@utu.fi

PAL can be installed by opening R and typing devtools::install_github("elolab/PAL") (requires package devtools to be installed). Notably, usage of PAL requires installation of R package PASI (devtools::install_github("elolab/PASI")) [2].

### Input and output

The only mandatory input from the user are **data** and **grouplabels**.
| Input | Description |
| ----------- | ----------- |
| data | A data frame of normalized gene expression data in log2 scale (rows: Entrez genes, cols: samples) |
| info | A data frame including all variables to be used in the analysis as columns, see details below. Rows correspond to samples in **data**. Rows and columns should be named. |
| grouplabels | Name of the column indicating sample groups in **info** |
| pathwayadress | A directory path to the user's own pathways (see format below) to be analysed on top of the KEGG pathways. If NULL (default), only the automatically accessed KEGG pathways are analysed. |
| useKEGG | TRUE (default) or FALSE indicating if KEGG pathways should be automatically accessed from API. Internet connection is required, if TRUE. |
| score | Either "activity" (default) or "deregulation" based on what the returned pathway scores should reflect (more details below). |
| nodemin | Indicates the minimum nuber of measured nodes in a pathway to be analysed (default 5). Pathways with fewer measured nodes are excluded from the analysis. |
| neutralize | NA (default) or name(s) of variable(s) (column(s) of **info**) that should be neutralised for prior to pathway analysis. The neutralised variables can be either numeric or character. |
| mainfeature | NA (default) or a name of a column in **info** containing either numeric values or character. If provided, pathways' significance levels according to this coefficient are returned. Can include NA for some samples.|
| seed | An integer used as a seed (default 1234) |

Data frame **info** should include minimum one column for **grouplabels**, which are indicated with integers (0 for controls). If the dataset contains no sample groups, use dummy value of 0 for all samples. Other columns can be either numeric or character. Notably, **info** can include also other columns than those indicated by **grouplabel**, **neutralize**, and **mainfeature**. Variables (columns) not indicated by any of those arguments are used in model fitting, but their effect is not neutralised. At this point redundant variables (e.g. gender does not change within a donor, so donor info can cover the variation caused bu gender) are not supported and cause singularity errors. Do not encode non-ordinal categorical variables with integers.

File format for the user defined pathways in **pathwayadress** is a .txt file and the directory should not include other .txt files than pathway files. Each pathway file should include the pathway's name as the first row and the following rows including nodes and relations. Node-lines should include only the Entrez gene id of the node. Relation-lines include three elemnts separated with a space: Entrez gene id of the start node (first), Entrez gene id of the end node (second), and either + or - indicating activation or inhibition, respectively (third). A toy example with five nodes (Entrez gene ids 5269, 8828, 10938, 1203, 8824) and three relations (5269 and 8828 activating 10938, and 8828 inhibiting 1203) is given below.

Pathway name here  
5269  
8828  
10938  
1203  
8824  
5269 10938 +  
8828 10938 +  
8828 1203 -

There are two options for argument **score**. The default one is "activity" and if it is selected, the final pathway scores indicate how active each pathway is in comparison to the other samples. Negative value indicate inactivity, value close to 0 normal activity, and positive value high activity. If **score** is set to "deregulation", the pathway scores indicate how normally the pathways behave as compared to a typical control sample. Value close to 0 means that the pathway is not deregulated and a high value means that it is.

Notably, this version of PAL can not neutralise for variables not overlapping between the sample groups defined in **grouplabels**. Therefore, for example donor can be used in model fitting, but can not be neutralised for.

The main function PAL returns a list with two or three elements: a matrix of pathway scores, a matrix of significance levels and coefficients (only if argument **mainfeature** is provided) of the analysed pathways, and an info matrix describing the analysed pathways. The pathway scores can be used for further analyses.

### Example

Run PAL for data frame mydata, which includes samples from two groups to obtain pathway activity scores.

	info = read.table("ExampleInfo.txt", sep="\t", quote="", stringsAsFactors=F)
	mydata = read.table("ExampleData.txt", sep="\t", quote="", stringsAsFactors=F)
	labels = rep(0, ncol(mydata))
	labels[grep("Case", colnames(mydata))] = 1
	timetodiagnosis = as.numeric(info[colnames(mydata), "TimeToDiagnosis"])
	res = PAL(mydata, grouplabels=labels, info=info[colnames(mydata),c("Age","Donor")], neutralize=c(T,F), mainfeature=timetodiagnosis)

##### References

[1] Add publication info when known

[2] Maria K Jaakkola, Aidan J McGlinchey, Riku Klén, and Laura L Elo: PASI: A novel pathway method to identify delicate group effects. PLoS one 2018; 13(7):e0199991. 	