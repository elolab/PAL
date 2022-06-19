\name{PAL}
\alias{PAL}

\title{
Pathway Analysis of Longitudinal data
}
\description{
Main function (and the only one for external use) of PAL performs sophisticated pathway analysis allowing complex study design, such as longitudinal data. The effect of some variables can be neutralised (see arguments 'info' and 'neutralise') prior to the pathway analysis.
}
\usage{
PAL(data, info, grouplabels, pathwayadress=NULL, useKEGG=TRUE, score="activity", nodemin=5, neutralize=NA, mainfeature=NA, seed=1234)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A data frame including gene expression data (rows: Entrez genes, cols: samples) as log2 values }
  \item{info}{ A data frame including all variables to be used in the analysis as columns, see details below. Rows correspond to samples in 'data'. Rows and columns should be named. }
  \item{grouplabels}{ Name of the column indicating sample groups in 'info' }
  \item{pathwayadress}{ A directory path to the user's own pathways (see format below) to be analysed on top of the KEGG pathways. If NULL (default), only the automatically accessed KEGG pathways are analysed.}
  \item{useKEGG}{ TRUE (default) or FALSE indicating if KEGG pathways should be automatically accessed from API. Internet connection is required, if TRUE.}
  \item{score}{ Either "activity" (default) or "deregulation" based on what the returned pathway scores should reflect}
  \item{nodemin}{ Pathways with less than 'nodemin' measured nodes will be excluded from the analysis}
  \item{neutralize}{ NA (default) or name(s) of variable(s) (column(s) of 'info') that should be neutralised for prior to pathway analysis. The neutralised variables can be either numeric or character. }
  \item{mainfeature}{ NA (default) or a name of a column in 'info' containing either numeric values or character. If provided, pathways' significance levels according to this coefficient are returned. Can include NA for some samples. }
  \item{seed}{ An integer used as a seed (default 1234)}
}
\details{ 
Data frame 'info' should include minimum one column for 'grouplabels', which are indicated with integers (0 for controls). If the dataset contains no sample groups, use dummy value of 0 for all samples. Other columns can be either numeric or character. Notably, 'info' can include also other columns than those indicated by '*grouplabel', 'neutralize', and 'mainfeature'. Variables (columns) not indicated by any of those arguments are used in model fitting, but their effect is not neutralised. At this point redundant variables (e.g. gender does not change within a donor, so donor info can cover the variation caused bu gender) are not supported and cause singularity errors. Do not encode non-ordinal categorical variables with integers.

File format for the user defined pathways in 'pathwayadress' is a .txt file and the directory should not include other .txt files than pathway files. Each pathway file should include the pathway's name as the first row and the following rows including nodes and relations. Node-lines should include only the Entrez gene id of the node. Relation-lines include three elemnts separated with a space: Entrez gene id of the start node (first), Entrez gene id of the end node (second), and either + or - indicating activation or inhibition, respectively (third). A toy example with five nodes (Entrez gene ids 5269, 8828, 10938, 1203, 8824) and three relations (5269 and 8828 activating 10938, and 8828 inhibiting 1203) is given below.

Pathway name here  
5269  
8828  
10938  
1203  
8824  
5269 10938 +  
8828 10938 +  
8828 1203 -

There are two options for argument 'score'. The default one is "activity" and if it is selected, the final pathway scores indicate how active each pathway is in comparison to the other samples. Negative value indicate inactivity, value close to 0 normal activity, and positive value high activity. If 'score' is set to "deregulation", the pathway scores indicate how normally the pathways behave as compared to a typical control sample. Value close to 0 means that the pathway is not deregulated and a high value means that it is.

Notably, this version of PAL can not neutralise for variables not overlapping between the sample groups defined in 'grouplabels'. Therefore, for example donor can be used in model fitting (i.e. appear as a column in 'info'), but can not be neutralised for (i.e. argument 'neutralize' can not include it).
}

\value{
Returns a list with two or three elements: 1) a matrix of pathway scores, 2) a matrix of significance levels and coefficients (only if argument 'mainfeature' is provided) of the analysed pathways, and 3) an info matrix describing the analysed pathways. The pathway scores can be used for further analyses.
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
# Run PAL for data frame exampleData, which includes samples from two groups to obtain pathway activity scores. Data frame exampleInfo contains columns 'Age', 'Donor', 'TimeToDiagnosis', and 'SampleGroup'.
# library(PAL)
# res = PAL(data=exampleData, info=exampleInfo, grouplabels="SampleGroup", neutralise="Age", mainfeature="TimeToDiagnosis")
}

\keyword{ design }
