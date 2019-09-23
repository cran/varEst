#' @title Variance Estimation with Ensemble method
#'
#' @description Estimation of error variance using ensemble method which combines bootstraping and sampling with srswor in ultrahigh dimensional dataset.
#'
#' @param x
#'
#' @param y
#'
#' @param a
#'
#' @param b
#'
#' @param d
#'
#' @return Error variance
#'
#' @examples
#'
#' @export
ensemble <- function(x,y,a,b,d,method= c("spam","lasso","lsr")){

  method = match.arg(method)


  p<- ncol(x)
  n<- nrow(x)


  if (method == "spam"){
    requireNamespace("SAM")
    spam_fit_n1 <- samQL(x,y,p=1)
    w_n1 <- row(as.matrix(spam_fit_n1$w[,30]))[which(spam_fit_n1$w[,30] != 0)]
    w_order_n1 <- head(order(spam_fit_n1$w[,30],decreasing = TRUE),d)
    w_value_n1 <- as.matrix(spam_fit_n1$w[,30])[w_order_n1,]
    spam_selected_feature_n1<- w_order_n1
    srswor_index <- combn(spam_selected_feature_n1,floor(d/2),FUN = NULL, simplify = FALSE)
    s <- dim(combn(d,floor(d/2)))[2]
    M1 <- dim(combn(d,floor(d/2)))[1]


    final_var <- vector("list",b)
    bs<- vector("list",b)
    xb <- vector("list",b)
    for (i in 1:b){
      bs[[i]] <- sample(1:n,n, replace = TRUE)
      xb[[i]] <- x[bs[[i]],]
      var <- numeric()
      selected_x <- vector("list",s)
      for (j in 1:s){
        selected_x[[j]] <- xb[[i]][,srswor_index[[j]]]

        fit_x2 <- lm(y ~ selected_x[[j]] - 1)
        var[[j]] <- sum((fit_x2$resid)^2)/(n - M1)
      }
      final_var[[i]] <- var
    }

    final_spam_var2 <- numeric()
    final_spam_var1<- vector("list",b)
    for (i in 1:b){
      final_spam_var<- numeric()
      for (j in 1:b){
        final_spam_var[[j]]<- mean(final_var[[i]][[j]])
      }
      final_spam_var1[[i]] <- final_spam_var
      final_spam_var2[[i]] <- mean(final_spam_var1[[i]])
    }

    final_spam_var3 <- mean(final_spam_var2)

    return(final_spam_var3)
  }


  if (method == "lasso"){
    requireNamespace("glmnet")
    lasso_fit_n1 <- glmnet(y=y,x=x,alpha= a,intercept=FALSE)
    lambda_n1=tail(lasso_fit_n1$lambda,1)
    beta_n1=coef(lasso_fit_n1,s=lambda_n1)
    beta_n1<-as.matrix(beta_n1)
    beta1_n1<- as.matrix(beta_n1[-1])
    select_beta_lasso_n1<- head(order(abs(beta1_n1), decreasing= TRUE),d)

    srswor_index <- combn(select_beta_lasso_n1,floor(d/2),FUN = NULL, simplify = FALSE)
    s <- dim(combn(d,floor(d/2)))[2]
    M1 <- dim(combn(d,floor(d/2)))[1]


    final_var <- vector("list",b)
    bs<- vector("list",b)
    xb <- vector("list",b)
    for (i in 1:b){
      bs[[i]] <- sample(1:n,n, replace = TRUE)
      xb[[i]] <- x[bs[[i]],]
      var <- numeric()
      selected_x <- vector("list",s)
      for (j in 1:s){
        selected_x[[j]] <- xb[[i]][,srswor_index[[j]]]

        fit_x2 <- lm(y ~ selected_x[[j]] - 1)
        var[[j]] <- sum((fit_x2$resid)^2)/(n - M1)
      }
      final_var[[i]] <- var
    }

    final_lasso_var2 <- numeric()
    final_lasso_var1<- vector("list",b)
    for (i in 1:b){
      final_lasso_var<- numeric()
      for (j in 1:b){
        final_lasso_var[[j]]<- mean(final_var[[i]][[j]])
      }
      final_lasso_var1[[i]] <- final_lasso_var
      final_lasso_var2[[i]] <- mean(final_lasso_var1[[i]])
    }
    final_lasso_var3 <- mean(final_lasso_var2)


    return(final_lasso_var3)
  }


  if (method == "lsr"){

    pValues_n<-numeric()
    for(i in 1:p){
      fitlm_n<- lm(y ~ x[,i])
      pValues_n[i]<-summary(fitlm_n)$coeff[2,4]
    }

    ranking_n<-order(pValues_n)
    pindex_n<- head(ranking_n,100)


    fitlm1_n<- lm(y~x[,pindex_n]-1)

    requireNamespace("lm.beta")
    lmbeta_n<- lm.beta(fitlm1_n)
    reg_beta_n<- lmbeta_n$standardized.coefficients
    reg_beta_n<- as.matrix(reg_beta_n)
    reg_beta_n[is.na(reg_beta_n)] <- 0
    reg_beta_select_n <- head(order(abs(reg_beta_n),decreasing = TRUE),d)
    reg_feature_name_n <- reg_beta_n[reg_beta_select_n,0]


    final_p_n<- summary(fitlm1_n)$coeff[,4]
    final_ranking_n<- order(final_p_n)
    final_p_n<- as.matrix(final_p_n)
    final_select_p_n <- row(final_p_n)[which(final_p_n<=0.05)]
    new_x <- x[,pindex_n]


    srswor_index <- combn(reg_beta_select_n,floor(d/2),FUN = NULL, simplify = FALSE)
    s <- dim(combn(d,floor(d/2)))[2]
    M1 <- dim(combn(d,floor(d/2)))[1]


    final_var <- vector("list",b)
    bs<- vector("list",b)
    xb <- vector("list",b)
    for (i in 1:b){
      bs[[i]] <- sample(1:n,n, replace = TRUE)
      xb[[i]] <- x[bs[[i]],]
      var <- numeric()
      selected_x <- vector("list",s)
      for (j in 1:s){
        selected_x[[j]] <- xb[[i]][,srswor_index[[j]]]

        fit_x2 <- lm(y ~ selected_x[[j]] - 1)
        var[[j]] <- sum((fit_x2$resid)^2)/(n - M1)
      }
      final_var[[i]] <- var
    }

    final_lsr_var2 <- numeric()
    final_lsr_var1<- vector("list",b)
    for (i in 1:b){
      final_lsr_var<- numeric()
      for (j in 1:b){
        final_lsr_var[[j]]<- mean(final_var[[i]][[j]])
      }
      final_lsr_var1[[i]] <- final_lsr_var
      final_lsr_var2[[i]] <- mean(final_lsr_var1[[i]])
    }

    final_lsr_var3 <- mean(final_lsr_var2)

    return(final_lsr_var3)
  }
}
