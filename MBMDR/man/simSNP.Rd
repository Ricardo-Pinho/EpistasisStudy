\name{simSNP}
\alias{simSNP}
\alias{simSNPcont}

\title{
  Simulated data of SNPs
}

\description{
  simSNP is a data.frame containing a bivariate outcome as in a case-control study, a continous covariate
  and 10 simulated SNPs. SNPs 1 and 2 have an interaction effect on the outcome. Penetrances 
  are specified as in model 1 of Ritchie et al.(2003) \cr

  simSNPcont is a data.frame containing a continuous normally distributed outcome and 10 simulated SNPs,
  where SNPs 1 and 2 have an interaction effect on the outcome.
}

\usage{
	data(simSNP)
	data(simSNPcont)
	}

\format{
  simSNP is a data frame with 400 observations on the following 12 variables.
  
  \tabular{ll}{
	  \code{Y}      \tab  Outcome: 1=case, 0=control  \cr
	  \code{X}      \tab  Age													\cr
    \code{snp.1}  \tab  SNP genotype data 1         \cr
    \code{snp.2}  \tab  SNP genotype data 2         \cr
    ...           \tab    .....                     \cr
    \code{snp.10} \tab  SNP genotype data 10        \cr
  }
  
  SNP values are coded as: 0=common homozygous, 1=heterozygous, 2=variant homozygous.
  	
  simSNPcont has a continous outcome and the same SNPs.
}

\references{
  Ritchie M.D., Hahn L.W., Moore J.H. (2003) \emph{Power of multifactor dimensionality reduction for detecting 
  gene-gene interactions in the presence of genotyping error, missing data, phenocopy, and genetic heterogeneity.}
  Genetic Epidemiology 24, 150-157
}
