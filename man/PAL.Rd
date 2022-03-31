\name{PAL}
\alias{PAL}

\title{
Pathway Analysis of Longitudinal data
}
\description{
Main function (and the only one for external use) of PAL performs sophisticated pathway analysis allowing complex study design, such as longitudinal data. The effect of some variables can be neutralised (see arguments **info** and **neutralize** prior to the pathway analysis.)
}
\usage{
PAL(data, grouplabels, pathwayadress=NULL, datatype="rnaseq", noisedefault="automatic", score="activity", nodemin=5, info=NA, neutralize=NA, mainfeature=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A data frame including gene expression data (rows: Entrez genes, cols: samples) as log2 values }
  \item{grouplabels}{ An integer vector indicating sample groups. Control samples should be indicated with label 0. }
  \item{pathwayadress}{ A directory path to the user's own pathways to be analysed on top of the KEGG pathways. If NULL (default), only the automatically accessed KEGG pathways are analysed.}
  \item{datatype}{ Either "microarray" or "rnaseq" (default) }
  \item{noisedefault}{ The default cutoff (numeric value or character ”automatic” (default)) for real signal. If set to NA, no lowly expressed genes are filtered out. }
  \item{score}{ Either "activity" (default) or "deregulation" based on what the returned pathway scores should reflect}
  \item{nodemin}{ Pathways with less than **nodemin** measured nodes will be excluded from the analysis}
  \item{info}{ A data frame including all variables to be used in the model fitting in the neutralisation step as columns. Rows correspond to samples in **data**. Rows and columns should be named. }
  \item{neutralize}{ A logical vector indicating which coefficients in **info** should be neutralised (T=neutralise, F=use in model fitting, but don't neutralise). The length should match the number of columns in **info**. }
  \item{mainfeature}{ A numeric or character vector corresponding to samples (cols in 'data'). If provided, pathways' significance levels according to this coefficient are returned. }
}
\details{ 
In case the expression data contains no sample groups, **grouplabels** can be set to dummy value of only zeros rep(0,ncol(data)).

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

Argument **datatype** defines the default value for **noisedefault** in log2 scale, which is different fot microarray data and RNAseq data (RNAseq has lower noise level). If the analysed data is neither, do nothing about argument **datatype** and set argument **noisedefault** to NA, and do the possible filtering of the data prior to PAL analysis.

If **noisedefault** is not NA, PAL filters out genes with median expression below a cutoff in all sample groups. The cutoff is detected from the data and **noisedefault** is the log2-scale default value for it.  

There are two options for argument **score**. The default one is "activity" and if it is selected, the final pathway scores indicate how active each pathway is in comparison to the other samples. Negative value indicate inactivity, value close to 0 normal activity, and positive value high activity. If **score** is set to "deregulation", the pathway scores indicate how normally the pathways behave as compared to a typical control sample. Value close to 0 means that the pathway is not deregulated and a high value means that it is.
}

\value{
Returns a list with two or three elements: 1) a matrix of pathway scores, 2) a matrix of significance levels and coefficients (only if argument **mainfeature** is provided) of the analysed pathways, and 3) an info matrix describing the analysed pathways. The pathway scores can be used for further analyses.
}

\references{
ADD REFERENCE
}

\author{
Maria K. Jaakkola
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
# Run PAL for data frame mydata, which includes samples from two groups to obtain pathway activity scores.
# info = read.table("ExampleInfo.txt", sep="\t", quote="", stringsAsFactors=F)
# mydata = read.table("ExampleData.txt", sep="\t", quote="", stringsAsFactors=F)
# labels = rep(0, ncol(mydata))
# labels[grep("Case", colnames(mydata))] = 1
# timetodiagnosis = as.numeric(info[colnames(mydata), "TimeToDiagnosis"])
# res = PAL(mydata, grouplabels=labels, info=info[colnames(mydata),c("Age","Donor")], neutralize=c(T,F), mainfeature=timetodiagnosis)
}

\keyword{ design }
