\name{ensemble}
\alias{ensemble}
\title{Variance Estimation with Ensemble method}
\description{
  Estimation of error variance using ensemble method which combines bootstraping and sampling with srswor in ultrahigh dimensional dataset.
}
\usage{
ensemble(x,y,a,b,d,method=c("spam","lasso","lsr"))
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
In this method, both bootstrapping and simple random sampling without replacement are combined to estimate error variance. Variables are selected using Sparse Additive Models (SpAM) or LASSO or least squared regression (lsr) from the original datasets and all possible samples of a particular size are taken from the selected variables set with simple random sampling without replacement. With these selected samples error variance is estimated from bootstrap samples of the original datasets using least squared regression method. Finally the average of all the estimated variances is considered as the final estimate of the error variance.
}
\value{
 \item{Error variance}{}
}
\references{Ravikumar, P., Lafferty, J., Liu, H. and Wasserman, L. (2009). Sparse additive models. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 71(5), 1009-1030
\cr Tibshirani, R. (1996). Regression shrinkage and selection via the Lasso. \emph{Journal of Royal Statistical Society}, 58, 267-288
}
\author{Sayanti Guha Majumdar <\email{sayanti23gm@gmail.com}>, Anil Rai, Dwijesh Chandra Mishra
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
var <- ensemble(marker,pheno,1,10,10,"spam")
}

\keyword{Ensemble}
