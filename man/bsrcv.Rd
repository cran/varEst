\name{bsrcv}
\alias{bsrcv}

\title{
Variance Estimation with Bootstrap-RCV
}
\description{
Estimation of error variance using Bootstrap-refitted cross validation method in ultrahigh dimensional dataset.
}
\usage{
bsrcv(x,y,a,b,d,method=c("spam","lasso","lsr"))
}

\arguments{
  \item{x}{a matrix of markers or explanatory variables, each column contains one marker and each row represents an individual.}
  \item{y}{a column vector of response variable.}
  \item{a}{value of alpha, range is 0<=a<=1 where, a=1 is LASSO penalty and a=0 is Ridge penalty.If variable selection method is LASSO then providing value to a is compulsory. For other methods a should be NULL.}
  \item{b}{number of bootstrap samples.}
  \item{d}{number of variables to be selected from x.}
  \item{method}{variable selection method, user can choose any method among "spam", "lasso", "lsr"}
}
\details{
In this method, bootstrap samples are taken from the original datasets and then RCV (Fan et al., 2012) method is applied to each of these bootstrap samples.
}
\value{\item{Error variance}{}
}
\references{
Fan, J., Guo, S., Hao, N. (2012).Variance estimation using refitted cross-validation in ultrahigh dimensional regression. \emph{Journal of the Royal Statistical Society}, 74(1), 37-65\cr Ravikumar, P., Lafferty, J., Liu, H. and Wasserman, L. (2009). Sparse additive models. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 71(5), 1009-1030\cr Tibshirani, R. (1996). Regression shrinkage and selection via the Lasso. \emph{Journal of Royal Statistical Society}, 58, 267-288
}
\author{
Sayanti Guha Majumdar <\email{sayanti23gm@gmail.com}>, Anil Rai, Dwijesh Chandra Mishra
}
\examples{
## data simulation
marker <- as.data.frame(matrix(NA, ncol =500, nrow = 200))
for(i in 1:500){
marker[i] <- sample(1:3, 200, replace = TRUE, prob = c(1, 2, 1))
}
pheno <- marker[,1]*1.41+marker[,2]*1.41+marker[,3]*1.41+marker[,4]*1.41+marker[,5]*1.41

pheno <- as.matrix(pheno)
marker<- as.matrix(marker)

## estimation of error variance
var <- bsrcv(marker,pheno,1,10,5,"lasso")
}

\keyword{ Bootstrap RCV }
