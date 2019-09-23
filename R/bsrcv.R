#' @title Variance Estimation with Bootstrap-RCV
#'
#' @description Estimation of error variance using Bootstrap-Refitted cross validation method in ultrahigh dimensional dataset.
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
bsrcv <- function(x,y,a,b,d,method= c("spam","lasso","lsr")){

  method = match.arg(method)


  p<- ncol(x)
  n<- nrow(x)


  bs<- vector("list",b)
  xb <- vector("list",b)
  for (i in 1:b){
    bs[[i]] <- sample(1:n,n, replace = TRUE)
    xb[[i]] <- x[bs[[i]],]
  }


  if (method == "spam"){
    spam_var_rcv <- function(x,y,d){
      p<- ncol(x)
      n<- nrow(x)

      k <- floor(n/2)
      x1 <- x[1:k, ]
      y1 <- y[1:k]
      x2 <- x[(k + 1):n, ]
      y2 <- y[(k + 1):n]
      n1 <- k
      n2 <- n-k




      requireNamespace("SAM")
      spam_fit_n1 <- samQL(x1,y1,p=1)
      w_n1 <- row(as.matrix(spam_fit_n1$w[,30]))[which(spam_fit_n1$w[,30] != 0)]
      w_order_n1 <- head(order(spam_fit_n1$w[,30],decreasing = TRUE),d)
      w_value_n1 <- as.matrix(spam_fit_n1$w[,30])[w_order_n1,]
      spam_selected_feature_n1<- w_order_n1

      M1 <- length(spam_selected_feature_n1)
      selected_x2 <- x2[,spam_selected_feature_n1]

      fit_x2 <- lm(y2 ~ selected_x2 - 1)
      var1 <- sum((fit_x2$resid)^2)/(n - k - M1)



      spam_fit_n2 <- samQL(x2,y2,p=1)
      w_n2 <- row(as.matrix(spam_fit_n2$w[,30]))[which(spam_fit_n2$w[,30] != 0)]
      w_order_n2 <- head(order(spam_fit_n2$w[,30],decreasing = TRUE),d)
      w_value_n2 <- as.matrix(spam_fit_n2$w[,30])[w_order_n2,]
      spam_selected_feature_n2<- w_order_n2

      M2 <- length(spam_selected_feature_n2)
      selected_x1 <- x1[,spam_selected_feature_n2]

      fit_x1 <- lm(y1 ~ selected_x1 - 1)
      var2 <- sum((fit_x1$resid)^2)/(k - M2)


      var_rcv <- (var1 + var2)/2
      return(var_rcv)
    }


    spam_var <- numeric()
    for (i in 1:b){
      spam_var[[i]] <- spam_var_rcv(xb[[i]],y,d)
    }
    return(mean(spam_var))
  }


  if (method == "lasso"){
    lasso_var_rcv <- function(x,y,a,d){
      p<- ncol(x)
      n<- nrow(x)

      k <- floor(n/2)
      x1 <- x[1:k, ]
      y1 <- y[1:k]
      x2 <- x[(k + 1):n, ]
      y2 <- y[(k + 1):n]
      n1 <- k
      n2 <- n-k




      requireNamespace("glmnet")
      lasso_fit_n1 <- glmnet(y=y1,x=x1,alpha= a,intercept=FALSE)
      lambda_n1=tail(lasso_fit_n1$lambda,1)
      beta_n1=coef(lasso_fit_n1,s=lambda_n1)
      beta_n1<-as.matrix(beta_n1)
      beta1_n1<- as.matrix(beta_n1[-1])
      select_beta_lasso_n1<- head(order(abs(beta1_n1), decreasing= TRUE),d)


      M1 <- length(select_beta_lasso_n1)
      selected_x2 <- x2[,select_beta_lasso_n1]

      fit_x2 <- lm(y2 ~ selected_x2 - 1)
      var1 <- sum((fit_x2$resid)^2)/(n - k - M1)


      lasso_fit_n2 <- glmnet(y=y2,x=x2,alpha= a,intercept=FALSE)
      lambda_n2=tail(lasso_fit_n2$lambda,1)
      beta_n2=coef(lasso_fit_n2,s=lambda_n2)
      beta_n2<-as.matrix(beta_n2)
      beta1_n2<- as.matrix(beta_n2[-1])
      select_beta_lasso_n2<- head(order(abs(beta1_n2), decreasing= TRUE),d)


      M2 <- length(select_beta_lasso_n2)
      selected_x1 <- x1[,select_beta_lasso_n2]

      fit_x1 <- lm(y1 ~ selected_x1 - 1)
      var2 <- sum((fit_x1$resid)^2)/(k - M2)


      var_rcv <- (var1 + var2)/2
      return(var_rcv)
    }


    lasso_var <- numeric()
    for (i in 1:b){
      lasso_var[[i]] <- lasso_var_rcv(xb[[i]],y,a,d)
    }

    return(mean(lasso_var))
  }


  if (method == "lsr"){
    lsr_var_rcv <- function(x,y,d){
      p<- ncol(x)
      n<- nrow(x)

      k <- floor(n/2)
      x1 <- x[1:k, ]
      y1 <- y[1:k]
      x2 <- x[(k + 1):n, ]
      y2 <- y[(k + 1):n]
      n1 <- k
      n2 <- n-k




      pValues_n1<-numeric()
      for(i in 1:p){
        fitlm_n1<- lm(y1~x1[,i])
        pValues_n1[i]<-summary(fitlm_n1)$coeff[2,4]
      }

      ranking_n1<-order(pValues_n1)
      pindex_n1<- head(ranking_n1,100)


      fitlm1_n1<- lm(y1~x1[,pindex_n1]-1)

      requireNamespace("lm.beta")
      lmbeta_n1<- lm.beta(fitlm1_n1)
      reg_beta_n1<- lmbeta_n1$standardized.coefficients
      reg_beta_n1<- as.matrix(reg_beta_n1)
      reg_beta_n1[is.na(reg_beta_n1)] <- 0
      reg_beta_select_n1 <- head(order(abs(reg_beta_n1),decreasing = TRUE),d)
      reg_feature_name_n1 <- reg_beta_n1[reg_beta_select_n1,0]


      final_p_n1<- summary(fitlm1_n1)$coeff[,4]
      final_ranking_n1<- order(final_p_n1)
      final_p_n1<- as.matrix(final_p_n1)
      final_select_p_n1 <- row(final_p_n1)[which(final_p_n1<=0.05)]
      new_x2 <- x2[,pindex_n1]


      M1 <- length(reg_beta_select_n1)
      selected_x2 <- new_x2[,reg_beta_select_n1]

      fit_x2 <- lm(y2 ~ selected_x2 - 1)
      var1 <- sum((fit_x2$resid)^2)/(n - k - M1)

      ## second set

      pValues_n2<-numeric()
      for(i in 1:p){
        fitlm_n2<- lm(y2~x2[,i])
        pValues_n2[i]<-summary(fitlm_n2)$coeff[2,4]
      }

      ranking_n2<-order(pValues_n2)
      pindex_n2<- head(ranking_n2,100)


      fitlm1_n2<- lm(y2~x2[,pindex_n2]-1)

      requireNamespace("lm.beta")
      lmbeta_n2<- lm.beta(fitlm1_n2)
      reg_beta_n2<- lmbeta_n2$standardized.coefficients
      reg_beta_n2<- as.matrix(reg_beta_n2)
      reg_beta_n2[is.na(reg_beta_n2)] <- 0
      reg_beta_select_n2 <- head(order(abs(reg_beta_n2),decreasing = TRUE),d)
      reg_feature_name_n2 <- reg_beta_n2[reg_beta_select_n2,0]


      final_p_n2<- summary(fitlm1_n2)$coeff[,4]
      final_ranking_n2<- order(final_p_n2)
      final_p_n2<- as.matrix(final_p_n2)
      final_select_p_n2 <- row(final_p_n2)[which(final_p_n2<=0.05)]
      new_x1 <- x1[,pindex_n2]


      M2 <- length(reg_beta_select_n2)
      selected_x1 <- new_x1[,reg_beta_select_n2]

      fit_x1 <- lm(y1 ~ selected_x1 - 1)
      var2 <- sum((fit_x1$resid)^2)/(n - k - M2)


      var_rcv <- (var1 + var2)/2
      return(var_rcv)

    }


    lsr_var <- numeric()
    for (i in 1:b){
      lsr_var[[i]] <- lsr_var_rcv(xb[[i]],y,d)
    }

    lsr_bs_var <- mean(lsr_var)
    return(lsr_bs_var)
  }
}
