#' @title Variance Estimation with kfold-RCV
#'
#' @description Estimation of error variance using k-fold Refitted cross validation in ultrahigh dimensional dataset.
#'
#' @param x
#'
#' @param y
#'
#' @param a
#'
#' @param k
#'
#' @param d
#'
#' @return Error variance
#'
#' @examples
#'
#' @export

 krcv <- function(x,y,a,k,d,method= c("spam","lasso","lsr")){

  method = match.arg(method)


  p<- ncol(x)
  n<- nrow(x)

  requireNamespace("caret")
  flds <- createFolds(y, k , list = TRUE, returnTrain = FALSE)
  split_up <- lapply(flds, function(ind, dat) dat[ind,], dat = y)
  for(i in 1:1010)
  {split_upx <- lapply(flds, function(i, dat) dat[i,], dat = x)}



  if (method == "spam"){

    w_n <- vector("list",k)
    w_order_n<- vector("list",k)
    w_value_n<- vector("list",k)
    spam_selected_feature_n<- vector("list",k)
    M<- numeric()

    requireNamespace("SAM")
    for (j in 1:k){
      spam_fit_n <- samQL(split_upx[[j]],split_up[[j]],p=1)

      w_n[[j]] <- row(as.matrix(spam_fit_n$w[,30]))[which(spam_fit_n$w[,30] != 0)]
      w_order_n[[j]] <- head(order(spam_fit_n$w[,30],decreasing = TRUE),d)
      w_value_n[[j]] <- as.matrix(spam_fit_n$w[,30])[w_order_n[[j]],]
      spam_selected_feature_n[[j]]<- w_order_n[[j]]
      M[[j]] <- length(spam_selected_feature_n[[j]])}

    final_selected_xM <- vector("list",k)
    final_var <- vector("list",k)
    for(t in 1:k){

      selected_xM <- vector("list",k)
      var <- numeric()

      for(l in 1:k){
        if (t != l){

          selected_xM[[l]] <- split_upx[[t]][,spam_selected_feature_n[[l]]]

          fit_xM <- lm(split_up[[t]] ~ selected_xM[[l]] - 1)
          var[[l]] <- sum((fit_xM$resid)^2)/(n/k - M[[t]])
        }
      }

      final_selected_xM[[t]] <- selected_xM
      final_var[[t]] <- var

    }

    var_krcv <- numeric()
    for(i in 1:k){
      var_krcv[[i]] <- mean(final_var[[i]], na.rm = TRUE)
    }

  }

  if (method == "lasso"){
    lambda_n <- numeric()
    beta_n<- vector("list",k)
    beta1_n<- vector("list",k)
    select_beta_lasso_n <- vector("list",k)
    M<- numeric()


    requireNamespace("glmnet")
    for (j in 1:k){
      lasso_fit_n <- glmnet(y=split_up[[j]],x=split_upx[[j]],alpha= a,intercept=FALSE)
      lambda_n <-tail(lasso_fit_n$lambda,1)
      beta_n<-coef(lasso_fit_n,s=lambda_n)
      beta_n<-as.matrix(beta_n)
      beta1_n[[j]]<- as.matrix(beta_n[-1])
      select_beta_lasso_n[[j]]<- head(order(abs(beta1_n[[j]]), decreasing= TRUE),d)


      M[[j]] <- length(select_beta_lasso_n[[j]])}

    final_selected_xM <- vector("list",k)
    final_var <- vector("list",k)
    for(t in 1:k){

      selected_xM <- vector("list",k)
      var <- numeric()

      for(l in 1:k){
        if (t != l){

          selected_xM[[l]] <- split_upx[[t]][,select_beta_lasso_n[[l]]]

          fit_xM <- lm(split_up[[t]] ~ selected_xM[[l]] - 1)
          var[[l]] <- sum((fit_xM$resid)^2)/(n/k - M[[t]])
        }
      }

      final_selected_xM[[t]] <- selected_xM
      final_var[[t]] <- var

    }

    var_krcv <- numeric()
    for(i in 1:k){
      var_krcv[[i]] <- mean(final_var[[i]], na.rm = TRUE)
    }

  }

  if (method == "lsr"){
    pindex_n <- vector("list",k)
    reg_beta_select_n <- vector("list",k)
    reg_feature_name_n <- vector("list",k)
    new_split_upx<- vector("list",k)
    M <- numeric()

    for (j in 1:k){
      pValues_n<-numeric()
      for(i in 1:p){
        fitlm_n<- lm(split_up[[j]] ~ split_upx[[j]][,i])
        pValues_n[i]<-summary(fitlm_n)$coeff[2,4]
      }

      ranking_n<-order(pValues_n)
      pindex_n[[j]]<- head(ranking_n,100)


      fitlm1_n<- lm(split_up[[j]]~split_upx[[j]][,pindex_n[[j]]]-1)

      requireNamespace("lm.beta")
      lmbeta_n<- lm.beta(fitlm1_n)
      reg_beta_n<- lmbeta_n$standardized.coefficients
      reg_beta_n<- as.matrix(reg_beta_n)
      reg_beta_n[is.na(reg_beta_n)] <- 0
      reg_beta_select_n[[j]] <- head(order(abs(reg_beta_n),decreasing = TRUE),d)
      reg_feature_name_n[[j]] <- reg_beta_n[reg_beta_select_n[[j]],0]


      final_p_n<- summary(fitlm1_n)$coeff[,4]
      final_ranking_n<- order(final_p_n)
      final_p_n<- as.matrix(final_p_n)
      final_select_p_n <- row(final_p_n)[which(final_p_n<=0.05)]
      new_split_upx[[j]] <- split_upx[[j]][,pindex_n[[j]]]


      M[[j]] <- length(reg_beta_select_n[[j]])}

    final_selected_xM <- vector("list",k)
    final_var <- vector("list",k)
    for(t in 1:k){

      selected_xM <- vector("list",k)
      var <- numeric()

      for(l in 1:k){
        if (t != l){

          selected_xM[[l]] <- new_split_upx[[t]][,reg_beta_select_n[[l]]]

          fit_xM <- lm(split_up[[t]] ~ selected_xM[[l]] - 1)
          var[[l]] <- sum((fit_xM$resid)^2)/(n/k - M[[t]])
        }
      }

      final_selected_xM[[t]] <- selected_xM
      final_var[[t]] <- var

    }

    var_krcv <- numeric()
    for(i in 1:k){
      var_krcv[[i]] <- mean(final_var[[i]], na.rm = TRUE)
    }

  }

  final_var_krcv <- mean(var_krcv)
  return(final_var_krcv)
}
