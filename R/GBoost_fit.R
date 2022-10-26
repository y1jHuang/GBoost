#' Group boosting fitting
#'
#' @param X matrix of predictors including any adjustment variables, which may not be grouped, of dimension n $/times$ p.
#' @param Y vector of the response variable, should be length n.
#' @param group vector containing an indicator of group membership of the p predictors in \code{X}
#' @param total_steps integer value of the maximum number of iterations to run the boosting algorithm
#' @param step_size a double value used to control how much to increase the estimate by at each step
#' @param adj_var integer number for the number of adjustment variables, 999 used to denote no adjustment variables
#' @param stop_tol double value used to determine what change in AIC/BIC/EBIC is small enough to stop iterating
#' @param gamma double value used to determine the penalty when using EBIC as stopping criteria
#' @param lasso_lambda double value of lambda used for lasso if regressing out the effect of adjustment variables
#' @param weighted default value is FALSE, set to TRUE if you want to use a weighted version to summarize within groups for the first stage
#' @return final_beta contains a vector with the final coefficient estimates, update1 is a vector that contains the coefficients from the first stage estimation, AIC1 contains the AIC path of the iterations from the first stage update, and BIC2 contains the BIC path from the second stage iterations
#' @export
GBoost_fit <- function(X, Y, group, total_steps=50000, step_size=c(1e-1, 1e-4), adj_var = 999, stop_tol=-1e-7, gamma = 0, lasso_lambda = 0.0314, weighted = 'n'){
  if(weighted=='n'){
    output = boosting_fitting_regress(X, Y, group, total_steps, step_size, adj_var, stop_tol, gamma, lasso_lambda)
    return(list(beta = output$final_beta, group_beta = output$update1, AIC = output$AIC1, BIC = output$BIC2, Z = output$Z))
  }else{
    output = boosting_fitting_weighted(X, Y, group, total_steps, step_size, adj_var, stop_tol, gamma, lasso_lambda, weighted)
    return(list(beta = output$final_beta, group_beta = output$update1, AIC = output$AIC1, BIC = output$BIC2, groups_weighted = output$groups.weighted))
  }
}

# weighted group boosting function
boosting_fitting_weighted <- function(X, Y, group, total_steps=50000, step_size=c(1e-1, 1e-4), adj_var = 999, stop_tol=-1e-7, gamma = 0, lasso_lambda = 0.0314, weighted='lasso'){
  num.groups <- length(unique(group))
  group.unique <- unique(group)
  sd.group <- NULL
  if(adj_var[1] != 999){
    fit.lasso1 = glmnet::glmnet(X, Y, penalty.factor=c(rep(0,length(adj_var)), rep(1, ncol(X) - length(adj_var))), lambda = lasso_lambda)
    fitted_values <- X[,adj_var]%*%fit.lasso1$beta[adj_var]
    res <- Y - fitted_values
    X.removed_adj = X[,-adj_var]  
    index.select <- which(!(c(1:ncol(X)) %in% adj_var))
  }else{
    X.removed_adj = X
    index.select <- c(1:ncol(X))
    res = Y
  }
  
  X.group <- matrix(NA, nrow = nrow(X.removed_adj), ncol = num.groups*2)
  sig <- c()
  group.vec.new <- rep(NA, ncol(X.removed_adj))
  group.size = rep(0, num.groups)
  sig.size = rep(0, 2*num.groups)
  for(g in 1:num.groups){
    group.tmp <- group.unique[g]
    idx.group <- which(group == group.tmp)
    # array of FC values of groups
    X.tmp <- X.removed_adj[,idx.group] # combine (sum or average) over groups
    group.size[g] = length(idx.group)
    if(is.vector(X.tmp)){
      X.group[,(2*g-1)] = X.tmp
      sig[2*g-1] = 1
      sig.size[2*g-1] = 1
    }else{
      if (weighted=='lasso') {
        mdl.tmp <- glmnet::cv.glmnet(X.tmp, res, alpha = 1) # fit lasso to the edges in the group
        beta.tmp <- as.vector(stats::coef(mdl.tmp, s=mdl.tmp$lambda.min))[-1]
      }else if(weighted=='boosting'){
        mdl.tmp <- model_fitting(X.tmp, res, x_size = rep(1, ncol(X.tmp)), 
                                 total_steps, step_size[1], adj_var, 
                                 stop_tol, gamma, stop_method = "AIC")
        beta.tmp <- mdl.tmp$beta
      }else{
        print('Please type in a correct weighted method: lasso / boosting')
        return()
      }
      # vectors of beta values of positive and negative signals
      idx.neg = which(beta.tmp < 0)
      sig.neg <- beta.tmp[idx.neg]
      X.neg = X.tmp[, idx.neg]
      idx.pos = which(beta.tmp > 0)
      sig.pos <- beta.tmp[idx.pos]
      X.pos = X.tmp[, idx.pos]
      group.vec.new[idx.group[which(beta.tmp == 0)]] <- paste(group.unique[g], ":none", sep = "")
      group.vec.new[idx.group[which(beta.tmp > 0)]] <- paste(group.unique[g], ":positive", sep = "")
      group.vec.new[idx.group[which(beta.tmp < 0)]] <- paste(group.unique[g], ":negative", sep = "")
      
      if(length(sig.pos) > 1){
        # value of strong positive signal
        sig.strongPos <- max(sig.pos)
        # index of strong positive signal
        idx.strongPos <- which.max(sig.pos)
        # index of weak positive signals
        idx.weakPos <- which(sig.pos<sig.strongPos)
        if(length(idx.weakPos) == 1){
          weighted.pos <- rowMeans(cbind(X.pos[,idx.weakPos], X.pos[,idx.strongPos]))
        }else{
          weighted.pos <- apply(X.pos[,idx.weakPos], 2, function(x){rowMeans(cbind(x, X.pos[,idx.strongPos]))})
        }
        X.group[,(2*g-1)] = rowMeans(cbind(weighted.pos, X.pos[,idx.strongPos]))
        sig[2*g-1] = sig.strongPos
        sig.size[2*g-1] = length(sig.pos)
      }else if(length(sig.pos) == 1){
        sig.strongPos <- max(sig.pos)
        idx.strongPos <- which.max(sig.pos)
        X.group[,(2*g-1)] <- X.pos
        sig[2*g-1] = sig.strongPos
        sig.size[2*g-1] = 1
      }else{
        X.group[,(2*g-1)] <- NA
        sig[2*g-1] <- NA
        sig.size[2*g-1] = 0
      }
      if(length(sig.neg) > 1){
        sig.strongNeg <- min(sig.neg) # quantile(sig.neg, 0.1)
        idx.strongNeg <- which.min(sig.neg) # index of strong negative signal
        idx.weakNeg <- which(sig.neg>sig.strongNeg)
        if(length(idx.weakNeg) == 1){
          weighted.neg <- rowMeans(cbind(X.neg[,idx.weakNeg], X.neg[,idx.strongNeg]))
        }else{
          weighted.neg <- apply(X.neg[,idx.weakNeg], 2, function(x){rowMeans(cbind(x, X.neg[,idx.strongNeg]))})
        }
        X.group[,2*g] = rowMeans(cbind(weighted.neg, X.neg[,idx.strongNeg]))
        sig[2*g] = sig.strongNeg
        sig.size[2*g] = length(sig.neg)
      }else if(length(sig.neg) == 1){
        sig.strongNeg <- min(sig.neg) # quantile(sig.neg, 0.1)
        idx.strongNeg <- which.min(sig.neg) # index of strong negative signal
        X.group[,2*g] <- X.neg
        sig[2*g] = sig.strongNeg
        sig.size[2*g] = 1
      }else{
        X.group[,2*g] <- NA
        sig[2*g] <- NA
        sig.size[2*g] = 0
      }
    }
  }

  X.std = matrix(NA, nrow = nrow(X.removed_adj), ncol = num.groups)
  group.unique.new <- rep(NA, num.groups)
  feat = c(":positive", ":negative", ":none")
  size = rep(0,num.groups)
  for (i in 1:num.groups) {
    X.tmp = X.group[, (2*i-1):(2*i)]
    sig.tmp = sig[(2*i-1):(2*i)]
    flag.tmp = is.na(sig.tmp)
    size.tmp = sig.size[(2*i-1):(2*i)]
    if (all(flag.tmp)) {
      group.unique.new[i] <- paste(group.unique[g], feat[3], sep = "")
      size[i] = 0
    }else{
      if (any(flag.tmp)) {
        X.std[,i] = X.tmp[, which(!flag.tmp)]
        group.unique.new[i] <- paste(group.unique[i], feat[which(!flag.tmp)], sep = "")
        size[i] = size.tmp[which(!flag.tmp)]
      }else{
        X.std[,i] = X.tmp[, which.max(abs(sig.tmp))]
        group.unique.new[i] <- paste(group.unique[i], feat[which.max(abs(sig.tmp))], sep = "")
        size[i] = size.tmp[which.max(abs(sig.tmp))]
      }
    }
  }
  idx.sig = !(apply(is.na(X.std), 2, any))
  X.shrink = X.std[,idx.sig]
  X.shrink <- apply(X.shrink, 2, scale, center = TRUE, scale = TRUE)
  time.start <- proc.time()
  boosting.group.fit <- model_fitting(X.shrink, scale(Y), x_size = rep(1, ncol(X.shrink)), total_steps, step_size[1], adj_var, stop_tol, gamma, stop_method = "AIC")
  time.groupsel = (proc.time() - time.start)[3]
  beta.pre = rep(0, num.groups)
  beta.pre[idx.sig] = boosting.group.fit$beta
  beta.update1 <- rep(0, length = ncol(X))
  for(g in 1:num.groups){
    beta.update1[index.select[which(group == group.unique[g])]] <- beta.pre[g]
  }
  if(adj_var[1] != 999){
    beta.update1[adj_var] <- 0
  }
  if(all(boosting.group.fit$beta == 0)){
    print("No groups selected.")
    X.red = X
    idx.keep = NULL
  }else{
    groups.selected <- group.unique[which(boosting.group.fit$beta != 0)]
    idx.keep <- which(group %in% groups.selected)
    if(adj_var[1] != 999){
      X.red <- as.matrix(X.removed_adj[,idx.keep])
    }else{
      X.red <- X.removed_adj[,idx.keep]
    }
  }
  step_factor = mean(group.size)
  boosting.ind.fit <- model_fitting(X.red, Y, x_size = rep(1, ncol(X.red)), total_steps*step_factor, step_size[2], adj_var, stop_tol, gamma, stop_method = "BIC")
  if(adj_var[1] != 999){
    beta.update2 <- rep(0, length = (ncol(X) - length(adj_var)))
    beta.update2[idx.keep] <- boosting.ind.fit$beta
  }else{
    beta.update2 <- rep(0, length = ncol(X))
    beta.update2[idx.keep] <- boosting.ind.fit$beta
  }
  return(list(final_beta = beta.update2, update1 = beta.update1, AIC1 = boosting.group.fit$AIC, BIC2 = boosting.ind.fit$BIC, groups.weighted = group.vec.new))
}
# function to fit group boosting algorithm and first regressing out the affect of adjustment variables with lasso
boosting_fitting_regress <- function(X, Y, group, total_steps=50000, step_size=c(1e-1, 1e-4), adj_var = 999, stop_tol=-1e-7, gamma = 0, lasso_lambda = 0.0314){
  num.groups <- length(unique(group))
  group.unique <- unique(group)
  group.size <- NULL
  sd.group <- NULL
  if(adj_var[1] != 999){
    fit.lasso1 = glmnet::glmnet(X, Y, penalty.factor=c(rep(0,length(adj_var)), rep(1, ncol(X) - length(adj_var))), lambda = lasso_lambda)
    fitted_values <- X[,adj_var]%*%fit.lasso1$beta[adj_var]
    residuals1 <- Y - fitted_values
    X.removed_adj = X[,-adj_var]
    index.select <- which(!(c(1:ncol(X)) %in% adj_var))
  }else{
    X.removed_adj = X
    index.select <- c(1:ncol(X))
    residuals1 = Y
  }
  X.group <- matrix(NA, nrow = nrow(X.removed_adj), ncol = num.groups)
  for(g in 1:num.groups){
    group.tmp <- group.unique[g]
    X.tmp <- X.removed_adj[,which(group == group.tmp)]
    if(is.vector(X.tmp)){
      X.group[,g] = X.tmp
    }else{
      X.group[,g] = rowMeans(X.tmp)
    }
    group.size = c(group.size, sum(group == group.unique[g]))
  }
  time.start <- proc.time()
  boosting.group.fit <- model_fitting(X.group, residuals1, x_size = rep(1, ncol(X.group)), total_steps, step_size[1], adj_var=adj_var, stop_tol, gamma, stop_method = "AIC")
  time.groupsel = (proc.time() - time.start)[3] # takes about 20 minutes for real data with 1000 steps
  if(adj_var[1] != 999){
    beta.update1 <- rep(0, length = (ncol(X) - length(adj_var)))
    for(g in 1:num.groups){
      beta.update1[which(group == group.unique[g])] <- boosting.group.fit$beta[g]
    }
  }else{
    beta.update1 <- rep(0, length = ncol(X))
    for(g in 1:num.groups){
      beta.update1[which(group == group.unique[g])] <- boosting.group.fit$beta[g]
    }
  }
  if(all(boosting.group.fit$beta == 0)){
    print("No groups selected.")
    X.red = X.removed_adj
  }else{
    groups.selected <- group.unique[which(boosting.group.fit$beta != 0)]
    idx.keep <- which(group %in% groups.selected)
    if(adj_var[1] != 999){
      X.red <- as.matrix(X.removed_adj[,idx.keep])
    }else{
      X.red <- X.removed_adj[,idx.keep]
    }
  }
  step_factor = mean(group.size)
  if(is.vector(X.red)){
    print("Only one variable selected")
    boosting.ind.fit <- model_fitting(X = X.red, y = residuals1, x_size = 1, total_steps = total_steps*step_factor, step_size[2], adj_var, stop_tol, gamma, stop_method = "BIC") # don't use residuals - use Y directly
  }else{
    boosting.ind.fit <- model_fitting(X = X.red, y = residuals1, x_size = rep(1, ncol(X.red)), total_steps = total_steps*step_factor, step_size[2], adj_var, stop_tol, gamma, stop_method = "BIC") # don't use residuals - use Y directly
  }
  if(adj_var[1] != 999){ # needed if adjustment variables
    beta.update2 <- rep(0, length = (ncol(X) - length(adj_var)))
    beta.update2[idx.keep] <- boosting.ind.fit$beta
  }else{
    beta.update2 <- rep(0, length = ncol(X))
    beta.update2[idx.keep] <- boosting.ind.fit$beta
  }
  return(list(final_beta = beta.update2, update1 = beta.update1, AIC1 = boosting.group.fit$AIC, BIC2 = boosting.ind.fit$BIC, Z = X.group))
}



# testing cpp function with example #### (put this in a different file)
# set.seed(2020)
# dat <- simul_dat_linear_model(n = 100, p = 10, rho_X = 0.5, effect_size = 0.5)
# step_size <- 0.01
# total_steps <- 5000
# stop_tol <- -1e-7
# beta_start <- rep(0, ncol(dat$X))
# adj_var = c(1,2,3)
# boosting_res = model_fitting(dat$X, dat$Y, rep(1,ncol(dat$X)), total_steps, step_size= 0.001, adj_var = 999, stop_tol, stop_method = "AIC")
# lasso_res = lasso_fitting(dat$X,dat$Y)
# enet_res = enet_fitting(dat$X,dat$Y,alpha=0.5)
# boosting_res$beta
# acctab = rbind(
#   compute_acc(boosting_res$beta,dat$beta),
#   compute_acc(lasso_res$beta,dat$beta),
#   compute_acc(enet_res$beta,dat$beta))
# rownames(acctab) = c("boosting no adj","lasso","enet")
# print(acctab,digit=2)
#
#
# compute_acc <- function(beta_est,true_beta){
#   true_idx = which(true_beta!=0)
#   select = factor(ifelse(beta_est!=0,1,0),levels=c(1,0))
#   true_select = factor(ifelse(true_beta!=0,1,0),levels=c(1,0))
#   tab = table(select,true_select)
#   TP = tab[1,1]
#   FP = tab[1,2]
#   FN = tab[2,1]
#   TN = tab[2,2]
#   Sensitivity  = TP/(TP+FN)
#   Specificity = TN/(TN+FP)
#   FDR = FP/(TP+FP)
#
#   return(c(null_mse = mean((beta_est[-true_idx]-true_beta[-true_idx])^2),
#            nonnull_mse = mean((beta_est[true_idx]-true_beta[true_idx])^2),
#            mse = mean((beta_est-true_beta)^2),
#            Sensitivity = Sensitivity,
#            Specificity = Specificity,
#            FDR = FDR))
# }

