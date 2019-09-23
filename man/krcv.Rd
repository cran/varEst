\name{krcv}
\alias{krcv}

\title{
Variance Estimation with kfold-RCV
}
\description{
Estimation of error variance using k-fold refitted cross validation in ultrahigh dimensional dataset.
}
\usage{
krcv(x,y,a,k,d,method=c("spam","lasso","lsr"))
}

\arguments{
  \item{x}{a matrix of markers or explanatory variables, each column contains one marker and each row represents an individual.}
  \item{y}{a column vector of response variable.}
  \item{a}{value of alpha, range is 0<=a<=1 where, a=1 is LASSO penalty and a=0 is Ridge penalty.If variable selection method is LASSO then providing value to a is compulsory. For other methods a should be NULL.}
  \item{k}{dataset is divided into this many numbers of sub-datasets.}
  \item{d}{number of variables to be selected from x.}
  \item{method}{variable selection method, user can choose any method among "spam", "lasso", "lsr"}
}
\details{
The error variance is estimated from a high dimensional datasets where number of parameters are more than number of individuals, i.e. p > n.k-fold RCV is an extended version of original RCV method (Fan et al.,  2012). In this case the datasets are divided into k equal size groups instead of 2 groups. Variables are selected using Sparse Additive Models (SpAM) or LASSO or least squared regression (lsr) from one group and variance is estimated using selected variables with ordinary least squared estimation from rest of the k-1 groups. Likewise, all the groups are covered and in the end, average value of all the variances from each group is the final error variance.
}
\value{
\item{Error variance}{}
}
\references{Fan, J., Guo, S., Hao, N. (2012).Variance estimation using refitted cross-validation in ultrahigh dimensional regression. \emph{Journal of the Royal Statistical Society}, 74(1), 37-65
\cr Ravikumar, P., Lafferty, J., Liu, H. and Wasserman, L. (2009). Sparse additive models. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 71(5), 1009-1030
\cr Tibshirani, R. (1996). Regression shrinkage and selection via the Lasso. \emph{Journal of Royal Statistical Society}, 58, 267-288
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
var <- krcv(marker,pheno,1,4,5,"spam")
}

\keyword{ kfold RCV }
