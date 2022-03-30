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

There are two options for argument **score**. The default one is "activity" and if it is selected, the final pathway scores indicate how active each pathway is in comparison to the other samples. Negative value indicate inactivity, value close to 0 normal activity, and positive value high activity. If **score** is set to "deregulation", the pathway scores indicate how normally the pathways behave as compared to a typical control sample. Value close to 0 means that the pathway is not deregulated and a high value means that it is.
.
}

\value{
Returns a list with two or three elements: 1) a matrix of pathway scores, 2) a matrix of significance levels and coefficients (only is argument **mainfeature** is provided) of the analysed pathways, and 3) an info matrix describing the analysed pathways. The pathway scores can be used for further analyses.
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
# labels = rep(c(0,1), times=c(20,18))
# timetodiagnosis = as.numeric(info[colnames(mydata), "TimeToDiagnosis"])
# res = PAL(mydata, grouplabels=labels, info=info[colnames(mydata),c("Age","Donor")], neutralize=c(T,F), mainfeature=timetodiagnosis)
}

\keyword{ design }
