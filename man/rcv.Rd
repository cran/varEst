\name{rcv}
\alias{rcv}
\title{Variance Estimation with Refitted Cross Validation(RCV)}
\description{
  Estimation of error variance using Refitted cross validation in ultrahigh dimensional dataset.
}
\usage{
rcv(x,y,a,d,method=c("spam","lasso","lsr"))
}
\arguments{
  \item{x}{a matrix of markers or explanatory variables, each column contains one marker and each row represents an individual.}
  \item{y}{a column vector of response variable.}
  \item{a}{value of alpha, range is 0<=a<=1 where, a=1 is LASSO penalty and a=0 is Ridge penalty. If variable selection method is LASSO then providing value to a is compulsory. For other methods a should be NULL.}
  \item{d}{number of variables to be selected from x.}
  \item{method}{variable selection method, user can choose any method among "spam", "lasso", "lsr"}
}
\details{The error variance is estimated from a high dimensional datasets where number of parameters are more than number of individuals, i.e. p > n. Refitted cross validation method (RCV) which is a two step method, is used to get the estimate of the error variance. In first step, dataset is divided into two sub-datasets and with the help of Sparse Additive Models (SpAM) or LASSO or least squared regression (lsr)  most significant markers(variables) are selected from the two sub-datasets. This results in two small sets of selected variables. Then using the set selected from 1st sub-dataset error variance is estimated from the 2nd sub-dataset with ordinary least square method and using the set selected from the 2nd sub-dataset error variance is estimated from the 1st sub-dataset with ordinary least square method. Finally the average of those two error variances are taken as the final estimator of error variance with RCV method.}
\value{\item{Error variance}{}
}
\references{Fan, J., Guo, S., Hao, N. (2012).Variance estimation using refitted cross-validation in ultrahigh dimensional regression. \emph{Journal of the Royal Statistical Society}, 74(1), 37-65\cr Ravikumar, P., Lafferty, J., Liu, H. and Wasserman, L. (2009). Sparse additive models. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 71(5), 1009-1030\cr Tibshirani, R. (1996). Regression shrinkage and selection via the Lasso. \emph{Journal of Royal Statistical Society}, 58, 267-288
}
\author{Sayanti Guha Majumdar <\email{sayanti23gm@gmail.com}>, Anil Rai, Dwijesh Chandra Mishra}
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
var <- rcv(marker,pheno,1,5,"spam")
}
\keyword{RCV}
