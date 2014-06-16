\name{mbmdr}
\Rdversion{2.1}
\alias{mbmdr}
\alias{print.mbmdr}

\title{
	Model Based Multifactor Dimensionality Reduction
}

\description{
  \code{mbmdr} implements the Model Based Multifactor Dimensionality Reduction (MB-MDR) method
  proposed by Calle et al.(2008) as a dimension reduction method for 
  exploring gene-gene interactions.
}
   
\usage{
  mbmdr(y, data, order, covar=NULL, exclude=NA, risk.threshold=0.1, 
        output=NULL, adjust=c("none","covariates","main effects","both"), 
        first.model=NULL, list.models=NULL, use.logistf=TRUE, 
        printStep1=FALSE, ...)
}

\arguments{

  \item{y}{
  	Vector containing the dependent variable.
  }

  \item{data}{
  	A data.frame (or object coercible by as.data.frame to a data frame) 
  	containing the SNP information with values 0,1,2. \cr 
  	For example: 0 = common homozygous genotype, 1 = heterozygous genotype, 
  	2 = variant homozygous genotype.
  }
 
  \item{order}{
  	Single integer that specifies the order of interactions to be analyzed.\cr
  	If \code{list.models = NULL} (value by default) all possible interactions of 
  	the specified order are analyzed.
  }
  
  \item{covar}{
  	\emph{(Optional)} A data.frame or object coercible by as.data.frame to a 
  	data frame containing the covariates for adjusting regression models. \cr 
  	Only used if \code{adjust="covariates"} or \code{adjust="both"}.
  }
 
  \item{exclude}{
  	\emph{(Optional)} Value/s of missing data. If missings in data are coded differently 
  	than \code{NA} it should be specified. 
  	For example \code{exclude=c(NA,-1)} specifies that both, \code{NA} and \code{-1}
  	indicate a missing value.
  }

  \item{risk.threshold}{
  	Threshold used at the first MB-MDR stage for defining the risk category of a multilocus
  	genotype. It should be a conservative value. 
  	The default value is \code{risk.threshold=0.1}.
  }

  \item{output}{
  	\emph{(Optional)} Output file name for storing \code{mbmdr} results on file, 
  	or \code{NULL} (default) for output as R object. If the number of models to be
    analyzed is too large, it is preferable to store the output in a file. This allows
    exploring partial results while mbmdr is still running and prevents from loosing all
    the information if R or \code{mbmdr} crashes during the process.
  }
 
  \item{adjust}{
  	Type of regressions adjustment. Options are \code{"none", "covariates", "main effects"}
  	or \code{"both"}. \cr
  	By default no adjustment is performed; \code{adjust="none"}.
  }
      	
  \item{first.model}{
  	\emph{(Optional)} Numerical vector of length equal to \code{order} 
    for specifying the initial interaction model for \code{mbmdr}; previous models will not be evaluated. 
    This is useful to continue \code{mbmdr} computation after a stop. \cr
    Note that, by default, \code{mbmdr} explores all possible interactions of a specified order. 
    If there are for example, 50 snps in data and \code{order=3}, \code{mbmdr} will start analyzing 
    the model (50,49,48), that means interaction between snps 50, 49 and 48 (column position 
    in data data.frame). The second model that \code{mbmdr} will analyze is (50,49,47). 
    After model (50,49,1), the next model will be (50,48,47), and the final model will be (3,2,1).\cr
 	  For example, if \code{mbmdr} stopped after the analysis of model (30,21,14), you can continue 
 	  the process by specifying \code{first.model=c(30,21,13)}.
	  Ids of snps must be in descended order.
    If \code{first.model=NULL} (by default) all models will be analyzed.
  }
   
  \item{list.models}{
  	\emph{(Optional)} Exhaustive list of models to be analyzed. Only models in list will be analyzed.
	  It can be: a vector of lenght \code{order} specifying a unique model, a matrix (\code{n x order}) 
	  containing models by rows, or a string for specifying a file with models by rows 
	  (all models must be of the same interaction order)\cr
    A \code{NULL} value (by default) indicates that all possible interactions will be analyzed.
  }

  \item{use.logistf}{
  	Boolean value indicating whether or not to use the \code{logistf} package 
    for logistic regressions when separation of data points is observed.\cr
    It only has effect if logistic regression (\code{family=binomial(link = "logit")}) is specified.
    (See \code{logistf} help for details). The default value is \code{TRUE}.
  }

  \item{printStep1}{
  	Boolean value. If true, the details of \code{mbmdr} step 1 are printed for every model.
    This slows the process, so it is only advisable when the number of models to analyze is small.\cr 
    By default \code{printStep1=FALSE}.
  }
  
  \item{...}{
  	For regression arguments: arguments to be passed to \code{glm} calls.\cr
  	Mainly to specify the error distribution and link function to be used in the regression models.\cr
  	For example, use \code{family=binomial(link=logit)} for specifying logistic regression or 
  	\code{gaussian(link = "identity")} for normal regression. \cr
  	(See \code{\link{family}} for details of family functions and \code{\link{glm}} for more options
  	of \code{glm} function).
  }

}

\details{
	MB-MDR is a method for identifying multi-locus genotypes that are associated with 
	a phenotype of interest, and allows to adjust for marginal and confounding effects.\cr

  The exploration of interacions is performed in three steps:\cr

  \emph{Step1}\cr
  Each genotype is tested for association with the response and is classified as high risk, low risk
  or not significant, and all genotypes of the same class are merged. The threshold for considering
  significant evidence is the value specified in \code{risk.threshold} (by default \code{risk.threshold=0.1}).\cr
  If \code{printStep1=TRUE}, the MBMDR function prints this classification.
  
  \emph{Step2}\cr
  For each risk categories, high and low, a new association test is performed. The result 
  provides a Wald statistic for the high and for the low categories.
  
  \emph{Step3}\cr
  The significance is explored through a permutation test on the maximum Wald statistics.

} 


\value{
  \code{mbmdr} returns an object of class \code{mbmdr} with the following attributes:
  
  \item{call}{
  The matched call.
  }
  \item{y}{
  The outcome used.
  }
  \item{data}{
  The SNPs data used.
  }
  \item{covar}{
  The covariate data used.
  }
  
  \item{result}{
  Dataframe with those interactions that have at least a significant genotype. 
  For each interaction (rows), the following information is returned:\cr
         
  \tabular{lll}{
    \code{SNP1...SNPx}\tab  \tab Names of snps in interaction.                                    \cr
    \code{NH   }      \tab  \tab Number of significant High risk genotypes in the interaction.    \cr
    \code{betaH}      \tab  \tab Regresion coeficient in step2 for High risk exposition.          \cr
    \code{WH   }      \tab  \tab Wald statistic for High risk category.                           \cr
    \code{PH   }      \tab  \tab P-value of the Wald test for the High risk category.             \cr
    \code{NL   }      \tab  \tab Number of significant Low risk genotypes in the interaction.     \cr
    \code{betaL}      \tab  \tab Regresion coeficient in step2 for Low risk exposition.           \cr
    \code{WL   }      \tab  \tab Wald statistic for Low risk category.                            \cr
    \code{PL   }      \tab  \tab P-value of the Wald test for the Low risk category.              \cr
    \code{MIN.P}      \tab  \tab Minimun p-value (\code{min(PH,PL)}) for the interaction model.   \cr
  }
  }

  If \code{printStep1} argument is set to \code{TRUE}, the result of the first step in 
  \code{mbmdr} is printed for each genotype with the following information:\cr
  
    \tabular{lll}{
      \code{...}      \tab  \tab Genotype.                                          \cr
      \code{cases}    \tab  \tab (only for case/control outcome) Number of cases with the specific genotype.        \cr
      \code{controls} \tab  \tab (only for case/control outcome) Number of controls with the specific genotype.     \cr
      \code{beta}     \tab  \tab Regression coefficient for this genotype.          \cr
      \code{p.value}     \tab  \tab Wald test p-value for this genotype.               \cr
      \code{category} \tab  \tab Predicted risk category for this genotype.         \cr
    }
}

\references{
  Calle M.L., Urrea V., Vellalta G., Malats N., Steen K.V. (2008) \emph{Improving strategies for 
  detecting genetic patterns of disease susceptibility in association studies.} 
  Statistics in Medicine 27, 6532-6546.
}


\author{
  Victor Urrea, Malu Calle, Kristel Van Steen, Nuria Malats
}


\examples{

#---  Dicotomous outcome  -------------

# load example data
data(simSNP)

# MB-MDR analysis of all possible 2nd order interactions (it may take some time)
# The order of the interactions to be explored is specified by order=2
# fit <- mbmdr(y=simSNP$Y,data=simSNP[,3:12],order=2,family=binomial(link=logit))
# print(fit)

# MB-MDR analysis of the epistatic effect of SNP1 and SNP2 (Model 2 1)
# The specific model to be analyzed is specified by list.models=c(2,1)
fit <- mbmdr(y=simSNP$Y,data=simSNP[,3:12],order=2,list.models=c(2,1),
             family=binomial(link=logit),printStep1=TRUE)
print(fit)


# MB-MDR analysis of the epistatic effect of SNP1 and SNP2, adjusted for variable X
# The specific model to be analyzed is specified by list.models=c(2,1)
# The adjustment statement is specified by adjust="covariates"
fit <- mbmdr(y=simSNP$Y,data=simSNP[,3:12],order=2,list.models=c(2,1),
             covar=simSNP$X,adjust="covariates",family=binomial(link=logit))
print(fit)


# MB-MDR analysis of all 2nd order interactions restricted to a subset of snps 
# (all possible 2nd order interactions among SNP1, SNP2 and SNP3 are explored,
# it may take some time)
# SNP1, SNP2 and SNP3 are placed in columns 3 to 5 of simSNP. This is specified
# by data=simSNP[,3:5]
# fit <- mbmdr(y=simSNP$Y,data=simSNP[,3:5],order=2,family=binomial(link=logit),
#             printStep1=TRUE)
# print(fit)


# MB-MDR analysis of all possible 3rd order interactions (it may take some time)
# The order of the interactions to be explored is specified by order=3
# fit <- mbmdr(y=simSNP$Y,data=simSNP[,3:12],order=3,family=binomial(link=logit))
# print(fit)


# MB-MDR analysis of the 3rd order epistatic effect of SNP1, SNP2 and SNP3
# The specific model to be analyzed is specified by list.models=c(3,2,1)
# fit <- mbmdr(y=simSNP$Y,data=simSNP[,3:12],order=3,list.models=c(3,2,1),
#              family=binomial(link=logit),printStep1=TRUE)
# print(fit)


#---  Continuous outcome  --------------
# load example data
data(simSNPcont)

# MB-MDR analysis of all possible 2nd order interactions (it may take some time)
# The order of the interactions to be explored is specified by order=2
# fit <- mbmdr(y=simSNPcont$Y,data=simSNPcont[,2:11],order=2)
# print(fit)

# MB-MDR analysis of the epistatic effect of SNP1 and SNP2 (Model 2 1)
# The specific model to be analyzed is specified by list.models=c(2,1)
fit <- mbmdr(y=simSNPcont$Y,data=simSNPcont[,2:11],order=2,
	           list.models=c(2,1),printStep1=TRUE)
print(fit)


}
