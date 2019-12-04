#' Estimate a partially functional linear regression model with latent group structures.
#' 
#' This is the proposed estimator, which is developed based on the functional principal component analysis.
#' The algorithm is a variant of the K-means clustering algorithm. 
#' The number of groups and the number of principal components are determined by the BIC criterion.
#' 
#' @param data_list A list of data. Several elements must be present in the list. The reponse \code{y},
#' the functional covariate \code{x_recv}, the scalar covariates \code{z}, and an index matrix \code{index}. 
#' The functional covariate \code{x_recv} must be generated from the \code{fda} package by, e.g., spline smoothing. 
#' The scalar covariates \code{z} is a matrix. The index matrix \code{index} is a data.frame recording the structure of the
#' data. The first column of \code{index} is the family number, the second column is the within family index. The column names of \code{index}
#' must be \code{ind_b} and \code{ind_w}.
#' 
#' @param num_group A vector of candidate number of groups.
#' 
#' @param num_pca A vector of candidate number of principal components.
#' 
#' @param est_fix_eff A logical value. If \code{TRUE}, then the fixed effects are estimated. Otherwise,
#' the fixed effects are not estimated
#' 
#' @param max_iter The maximum number of iteration. Default to be 100.
#' 
#' @param loc_search A logical value. If \code{TRUE}, conduct a local search to decrease the objective function.
#' 
#' @param group_index An initial value of the group membership, optional.
#' 
#' 
#' 
#' 
#' @export
bic_kmean_est <- function(data_list, num_group = 2, num_pca = 5,
                          est_fix_eff=TRUE, max_iter = 100, loc_search = FALSE, group_index = NULL){
  para_df  = expand.grid(num_group, num_pca)
  n        = nrow(data_list$index) # total number of observations
  n_b   = length(unique(data_list$index[,1]))
  #saveRDS(data_list, 'err_data.rds')
  est_treat = FALSE


    bic_factor  = 1
  

  # store the estimators
  est1 = c();    bicValue1 = 1e20;  lambdabic1 = c();
  
  a_vector    = c()
  b_vector    = c()
  sumq_vector = c()

  for( iter in 1:nrow(para_df)){
    num_group_t = para_df[iter, 1]
    num_pca_t   = para_df[iter, 2]
    #print(sprintf('num_group %d num_pca %d', num_group_t, num_pca_t))
    #cv_temp    = kfold_cv(data_list, K=K, num_group = num_group_t, num_pca = num_pca_t, nbasis = nbasis, gamma_basis=gamma_basis)

    if(loc_search){
      cv_temp1    = fdagroup_est_scale(data_list, num_group = num_group_t, num_pca = num_pca_t, est_fix_eff=est_fix_eff, est_treat = est_treat,
                                        max_iter = max_iter, scale_kmean = scale_kmean)
      cv_temp     = local_search(data_list, cv_temp1, num_group = num_group_t, num_pca = num_pca_t, est_fix_eff=est_fix_eff, est_treat = est_treat,
                                 max_iter = 10, scale_kmean = scale_kmean)
    } else {
      cv_temp    = fdagroup_est_scale(data_list, num_group = num_group_t, num_pca = num_pca_t, est_fix_eff=est_fix_eff, est_treat = est_treat,
                                      max_iter = max_iter, scale_kmean = FALSE, group_index = group_index)
    }
    # cat("iteration steps:")
    # cat(cv_temp$iter)
    # cat('\n')
    sumq        = ncol(data_list$z) + num_group_t * num_pca_t + as.numeric(est_treat) + as.numeric(est_fix_eff) * n_b

    tbicValue   = cv_temp$aa + bic_factor *log(n) * (sumq)/n  #*log(log(n+p))
    tbicValue[1] =   cv_temp$aa + 2 * (sumq)/n               # AIC

    a_vector    = c(a_vector, cv_temp$aa)
    b_vector    = c(b_vector, log(n) * (sumq)/n)
    sumq_vector = c(sumq_vector, sumq)

    for(i in 1:length(bic_factor)){
      if(tbicValue[i] < eval(parse(text = paste("bicValue", i, sep = "")))){
        assign(paste("est", i, sep = ""), cv_temp)
        assign(paste("bicValue", i, sep = ""), tbicValue[i])
        assign(paste("lambdabic", i, sep = ""), c(num_pca_t, num_group_t))
      }
    }
  }

  #cv_temp     = local_search(data_list, est8, num_group = lambdabic8[2], num_pca = lambdabic8[1], nbasis = nbasis,
  #                           gamma_basis=gamma_basis, est_fix_eff=est_fix_eff, est_treat = est_treat,
  #                           max_iter = 5, scale_kmean = scale_kmean)
  #est_8 = cv_temp

  # the parameter combination which minimize the cv-error.
  return(list(est1 = est1, bicValue1 = bicValue1, lambdabic1 = lambdabic1))
}

empty_group <-  function(group_index_old, num_group = NULL){
  # randomly assign subjects into empty groups.
  group_index = group_index_old
  uni_group   = unique(group_index)
  if(length(uni_group) < num_group){
    for(i in 1:10){ # try several times, incase some groups are empty again after the reassignment
      sub_assign = sample(1:length(group_index), num_group - length(uni_group), replace = FALSE)
      for(sub in 1:length(sub_assign)){
        group_index[sub_assign[sub]] = length(uni_group) + sub
      }
      if(length(unique(group_index)) == num_group){ # succeeded
        group_index = adjust_group(group_index)
        return(group_index)
      } else { # failed, try to assign again
        group_index = group_index_old
      }
    }
    return(group_index)
  } else {
    return(group_index_old)
  }
}

# used to expand x_pca_socre, such that it has columns ncol(x_pca_score)*length(group_index)
x_pca_m     = function(x_pca_score, group_index, index){

  # if only one group, return directly
  if(length(unique(group_index)) == 1){
    return(x_pca_score)
  }

  res = matrix(0, nrow(x_pca_score), ncol(x_pca_score)*length(unique(group_index)))
  n_b = length(unique(index[,1])) # group_index is of length n_b
  if(n_b != length(group_index)){
    stop('number of subjects is not the same as group_index')
  }
  group_df = data.frame(unique(index[,1]), group_index)
  colnames(group_df) = c('ind_b', 'group')
  group_long = left_join(index, group_df, by='ind_b')
  group_long$group = factor(group_long$group)

  # using formulas and model.matrix to replace the inefficient code
  colnames(x_pca_score) = paste('X', 1:ncol(x_pca_score), sep='')
  fmla   = as.formula(paste('~', paste("group:", paste('X', 1:ncol(x_pca_score), sep=''), sep = '', collapse = '+'), '+0', sep = ''))
  data_g = data.frame(group_long, x_pca_score)
  res_t  = model.matrix(fmla, data = data_g)

  # sort according to group names
  res_name = sort(colnames(res_t))
  res      = res_t[, res_name]
  colnames(res) = NULL

  # g_ind_uni_sort =  dplyr::dense_rank(sort(unique(group_index)))
  # for(i in 1:length(g_ind_uni_sort)){
  #   res[,((i-1)*ncol(x_pca_score)+1):((i)*ncol(x_pca_score))] = diag(group_long$group == g_ind_uni_sort[i]) %*% x_pca_score
  # }
  return(res)
}


#' @export
fdagroup_est_scale <- function(data_list, num_group = 2, num_pca = 5, est_fix_eff=TRUE,
                               est_treat=TRUE, max_iter = 100, group_index = NULL, subj_center = FALSE, scale_kmean = FALSE){
  # data. x should be discrete observation on the grid points t_grid. nbasis, the number of basis in smooth x.
  y      = data_list$y
  #x_dis  = data_list$x_dis
  x_recv = data_list$x_recv
  #t_grid = data_list$t_grid
  z      = data_list$z
  index  = data_list$index
  n_b    = length(unique(index[,1]))
  
  # FPCA
  x_recv_pca  = fda::pca.fd(x_recv, nharm = num_pca)
  x_pca_score = x_recv_pca$scores

  # If we center the observations subject-wise, then we do not have to estimate the ground mena mu and the fixed effects
  if(subj_center){
    temp_dat = data.frame(y, z, x_pca_score, index[,3], ind = index[,1]) # index[,1] has the subject information.
    temp_dat = as.data.frame(temp_dat %>% group_by(ind) %>% mutate_all(.funs = scale, scale = FALSE))
    y = temp_dat[,1]
    z = as.matrix(temp_dat[,2:(1+ncol(z))])
    x_pca_score = as.matrix(temp_dat[, (2+ncol(z)):(1+ ncol(z) + ncol(x_pca_score))])
    index[,3]   = temp_dat[, (2+ ncol(z) + ncol(x_pca_score))]
  }

  # scale the x_pca_score.
  if(scale_kmean){
    scale_temp = scale(x_pca_score, center = FALSE)
    scale_vec  = attr(scale_temp, 'scaled:scale')
    x_pca_score_scale = scale_temp
  } else {
    x_pca_score_scale = x_pca_score
    scale_vec         = rep(1, ncol(x_pca_score))
  }



  # a subject-wise summary matrix
  prodKron   = matrix(0,nrow=nrow(index),ncol=n_b)
  index_iter = 0
  struT      = (index %>% dplyr::group_by(ind_b) %>% dplyr::summarise(n_w=n()))$n_w
  for(i in 1:n_b ){
    prodKron[(index_iter+1):(index_iter+struT[i]),i] = 1;
    index_iter = index_iter + struT[i];
  }

  # iterating k-means.
  # initial value for the group index
  if(is.null(group_index)){
    lmfitz = lm(y~z)
    epi_df = data.frame(index, res=lmfitz$residuals^2) %>% dplyr::group_by(ind_b) %>% dplyr::summarize(resi2_mn=mean(res))
    group_index =  adjust_group(dplyr::dense_rank(kmeans(epi_df$resi2_mn, centers = num_group, nstart = 100)$cluster))
  }

  #max_iter = 100
  group_index_old = group_index
  zhat            = rep(0, ncol(z))
  xhat_subject    = matrix(0, num_pca, n_b)
  xhat_subject_old= xhat_subject
  xhat_group      = matrix(0, num_pca, num_group) # store the coefficient for each group,
  xhat_group_old  = xhat_group
  des_dummy       = create_design(index, est_fix_eff=est_fix_eff, mu = TRUE, est_treat = est_treat) # pseudo design, include fixed effects and treatment effects
  dummyhat        = rep(0, ncol(des_dummy))  # +1 to include ground mean
  # number of columns (groups) could less than num_group
  for(iter in 1:max_iter){
    # step 1 update group index
    if(iter >1){
      group_index_old = group_index
      resi            = (as.vector(y) - as.vector(des_dummy %*% dummyhat) - as.numeric(as.matrix(z %*% zhat)) - x_pca_score_scale %*% xhat_group)
      group_index     = apply(t(prodKron) %*% (resi^2), 1, which.min)
      group_index     = adjust_group(dplyr::dense_rank(group_index))

      # randomly assign subjects into empty groups. The returned groups are adjusted.
      group_index     = empty_group(group_index, num_group)
    }

    # step 2. estimation
    x_score_group   = x_pca_m(x_pca_score_scale, group_index, index) # for functional covariates
    X               = cbind(des_dummy, z, x_score_group)    # pseduo design matrix
    lmfit           = lm(y ~ X + 0)       # the intercept is added automatically.
    coef_fit        = coef(lmfit)
    coef_fit[which(is.na(coef_fit))] = 0
    dummyhat        = coef_fit[1:(ncol(des_dummy))]
    zhat            = coef_fit[(length(dummyhat)+1):(ncol(z) + length(dummyhat))]
    xhat_group_old  = xhat_group
    xhat_group      = diag(1/scale_vec) %*% matrix(coef_fit[(ncol(z)+length(dummyhat)+1):length(coef_fit)], nrow=num_pca)

    #
    uni_group_index  = (unique(group_index))
    xhat_subject_old = xhat_subject
    for(ix in 1:n_b){
      xhat_subject[,ix] = xhat_group[, which(group_index[ix] == uni_group_index)]
    }

    # stop iteration
    if( (iter>1 & identical(group_index, group_index_old)) | base::norm(xhat_subject-xhat_subject_old, type = 'F')<1e-5){
      break
    }
  }

  # value of the object function
  aa = log(sum(lmfit$residuals^2)/length(y))

  # represent estimation functional data form
  xhat_group_fd = fd(coef=x_recv_pca$harmonics$coefs %*% xhat_group, basisobj = x_recv_pca$harmonics$basis)
  xhat_sub_fd   = fd(coef=x_recv_pca$harmonics$coefs %*% xhat_subject, basisobj = x_recv_pca$harmonics$basis)

  # gound mean, fixed effects and treatment effects
  muhat = as.numeric(dummyhat[1])
  if(est_fix_eff){
    fixeffhat = c(0,as.numeric(dummyhat[2:n_b])) # now the fixed effects has length n_b, and the first element is 0
  } else {
    fixeffhat = rep(0, n_b)
  }
  if(est_treat){
    treathat = as.numeric(dummyhat[length(dummyhat)])
  } else {
    treathat = 0
  }
  return(list(muhat = muhat, fixeffhat = fixeffhat, treathat = treathat, zhat=zhat,
              group=group_index, xhat_subject = xhat_subject, xhat_sub_fd=xhat_sub_fd, aa=aa, iter = iter))
}


fdagroup_est_search <- function(data_list, num_group = 2, num_pca = 5, nbasis = 20, gamma_basis=NULL, est_fix_eff=TRUE,
                                est_treat=TRUE, max_iter = 100, group_index = NULL, subj_center = FALSE, scale_kmean = FALSE,
                                num_try = 100){
  # data. x should be discrete observation on the grid points t_grid. nbasis, the number of basis in smooth x.
  y      = data_list$y
  #x_dis  = data_list$x_dis
  x_recv = data_list$x_recv
  #t_grid = data_list$t_grid
  z      = data_list$z
  index  = data_list$index
  n_b    = length(unique(index[,1]))

  # smooth x
  if(is.null(gamma_basis)){
    gamma_basis = fda::create.fourier.basis(nbasis = 2*nbasis+1, dropind = c(1, seq(2, 2*nbasis, by=2)))
  }
  #x_recv  = fda::smooth.basis(argvals = t_grid, y = x_dis, fdParobj = gamma_basis)

  # FPCA
  x_recv_pca  = fda::pca.fd(x_recv, nharm = num_pca)
  x_pca_score = x_recv_pca$scores

  # scale the x_pca_score.
  if(scale_kmean){
    scale_temp = scale(x_pca_score, center = FALSE)
    scale_vec  = attr(scale_temp, 'scaled:scale')
    x_s = scale_temp
  } else {
    x_s = x_pca_score
    scale_vec = rep(1, ncol(x_pca_score))
  }

  # If we center the observations subject-wise, then we do not have to estimate the ground mena mu and the fixed effects
  if(est_fix_eff){
    temp_dat = data.frame(y, z, x_s, index[,3], ind = index[,1]) # index[,1] has the subject information.
    temp_dat = as.data.frame(temp_dat %>% group_by(ind) %>% mutate_all(.funs = scale, scale = FALSE))
    y_c = temp_dat[,1]
    z_c = as.matrix(temp_dat[,2:(1+ncol(z))])
    x_s_c = as.matrix(temp_dat[, (2+ncol(z)):(1+ ncol(z) + ncol(x_s))])
    treat_ind_c = temp_dat[, (2+ ncol(z) + ncol(x_s))]
  } else {
    y_c = y
    z_c = z
    x_s_c = x_s
    treat_ind_c = index[,3]
  }

  # a subject-wise summary matrix
  prodKron   = matrix(0,nrow=nrow(index),ncol=n_b)
  index_iter = 0
  struT      = (index %>% dplyr::group_by(ind_b) %>% dplyr::summarise(n_w=n()))$n_w
  for(i in 1:n_b ){
    prodKron[(index_iter+1):(index_iter+struT[i]),i] = 1;
    index_iter = index_iter + struT[i];
  }

  prodKron_mean   = prodKron %*% diag(1/struT)


  if(est_treat){
    z_treat_c = cbind(treat_ind_c, z_c)
    z_treat   = cbind(index[,3], z)
  } else {
    z_treat_c = z_c
    z_treat   = z
  }

  z_treat_hat     = rep(0, ncol(z_treat))
  xhat_subject    = matrix(0, num_pca, n_b)
  xhat_subject_old= xhat_subject
  xhat_group      = matrix(0, num_pca, num_group) # store the coefficient for each group,
  xhat_group_old  = xhat_group
  des_dummy       = create_design(index, est_fix_eff=est_fix_eff, mu = FALSE, est_treat = FALSE) # pseudo design, include fixed effects and treatment effects
  if(is.null(des_dummy)){
    des_dummy = create_design(index, est_fix_eff=TRUE, mu = FALSE, est_treat = FALSE)
    dummyhat  = rep(0, n_b)
  } else {
    dummyhat  = rep(0, ncol(des_dummy))  # the fixed effects only
  }

  z_treat_hat_best  = z_treat_hat
  xhat_subject_best = xhat_subject
  dummyhat_best     = dummyhat
  group_index_best  = rep(1, n_b)
  obj_val_best    = 1e100

  # iterating k-means.
  # initial value for the group index
  group_index_matrix = c()
  if(is.null(group_index)){
    lmfitz = lm(y_c ~ z_treat_c)
    epi_df = data.frame(index, res=lmfitz$residuals^2) %>% dplyr::group_by(ind_b) %>% dplyr::summarize(resi2_mn=mean(res))
    group_index =  adjust_group(dplyr::dense_rank(kmeans(epi_df$resi2_mn, centers = num_group, nstart = 100)$cluster))
  }
  if(num_group > 1){
    group_index_matrix = rbind(group_index, matrix(sample(1:num_group, size = n_b * num_try, replace = TRUE), nrow = num_try, ncol = n_b))
  } else {
    group_index_matrix = matrix(group_index, nrow=1)
  }

  for(iter_g in 1:nrow(group_index_matrix)){
    group_index     = adjust_group(group_index_matrix[iter_g,])
    group_index_old = group_index

    for(iter in 1:max_iter){
      # step 1 update group index
      if(iter >1){
        group_index_old = group_index
        resi            = (as.vector(y_c) - as.vector(des_dummy %*% dummyhat) - as.numeric(as.matrix(z_treat_c %*% z_treat_hat)) - x_s_c %*% xhat_group)
        group_index     = apply(t(prodKron) %*% (resi^2), 1, which.min)
        group_index     = adjust_group(dplyr::dense_rank(group_index))

        # randomly assign subjects into empty groups. The returned groups are adjusted.
        group_index     = empty_group(group_index, num_group)
      }

      # step 2. estimation
      x_score_group   = x_pca_m(x_s_c, group_index, index) # for functional covariates
      X               = cbind(z_treat_c, x_score_group)    # pseduo design matrix
      lmfit           = lm(y_c ~ X + 0)       # the intercept is added automatically.
      coef_fit        = coef(lmfit)
      coef_fit[which(is.na(coef_fit))] = 0
      z_treat_hat     = coef_fit[1:ncol(z_treat_c)]
      xhat_group_old  = xhat_group
      xhat_group      = diag(1/scale_vec) %*% matrix(coef_fit[(ncol(z_treat_c)+1):length(coef_fit)], nrow=num_pca)
      if(est_fix_eff){
        dummyhat        = c(t(prodKron_mean) %*% as.vector(y) - (t(prodKron_mean) %*% x_pca_m(x_pca_score, group_index, index) ) %*% c(xhat_group) - (t(prodKron_mean) %*% z_treat) %*% z_treat_hat)   # the fixed effects estimate
      } else {
        dummyhat = rep(0, n_b)
      }
      #
      uni_group_index  = (unique(group_index))
      xhat_subject_old = xhat_subject
      for(ix in 1:n_b){
        xhat_subject[,ix] = xhat_group[, which(group_index[ix] == uni_group_index)]
      }

      # stop iteration
      if( (iter>1 & identical(group_index, group_index_old)) | base::norm(xhat_subject-xhat_subject_old, type = 'F')<1e-5){
        break
      }
    }
    # the value of the object function
    obj_val = sum((as.vector(y) - as.vector(des_dummy %*% dummyhat) - as.numeric(as.matrix(z_treat %*% z_treat_hat)) - x_pca_m(x_pca_score, group_index, index) %*% c(xhat_group))^2)
    if(obj_val < obj_val_best){
      z_treat_hat_best  = z_treat_hat
      xhat_subject_best = xhat_subject
      dummyhat_best     = dummyhat
      group_index_best  = group_index
      obj_val_best      = obj_val
    }
  }

  # value of the object function
  #aa = log(sum(lmfit$residuals^2)/length(y))
  aa  = obj_val_best / length(y)

  # represent estimation functional data form
  #xhat_group_fd = fd(coef=x_recv_pca$harmonics$coefs %*% xhat_group, basisobj = x_recv_pca$harmonics$basis)
  xhat_sub_fd   = fd(coef=x_recv_pca$harmonics$coefs %*% xhat_subject_best, basisobj = x_recv_pca$harmonics$basis)

  # gound mean, fixed effects and treatment effects
  muhat = 0
  if(est_fix_eff){
    fixeffhat = dummyhat_best
  } else {
    fixeffhat = rep(0, n_b)
  }
  if(est_treat){
    treathat = z_treat_hat_best[1]
    zhat     = z_treat_hat_best[2:length(z_treat_hat_best)]
  } else {
    treathat = 0
    zhat     = z_treat_hat
  }
  return(list(muhat = muhat, fixeffhat = fixeffhat, treathat = treathat, zhat=zhat,
              group=group_index_best, xhat_subject = xhat_subject_best, xhat_sub_fd=xhat_sub_fd, aa=aa, iter = iter))
}


local_search <- function(data_list, est_res, num_group = 2, num_pca = 5, nbasis = 20, gamma_basis=NULL, est_fix_eff=TRUE,
                         est_treat=TRUE, max_iter = 10, group_index = NULL, subj_center = FALSE, scale_kmean = FALSE){
  # based on the return of the last function.
  y      = data_list$y
  #x_dis  = data_list$x_dis
  x_recv = data_list$x_recv
  #t_grid = data_list$t_grid
  z      = data_list$z
  index  = data_list$index
  n_b    = length(unique(index[,1]))

  # smooth x
  if(is.null(gamma_basis)){
    gamma_basis = fda::create.fourier.basis(nbasis = 2*nbasis+1, dropind = c(1, seq(2, 2*nbasis, by=2)))
  }
  #x_recv  = fda::smooth.basis(argvals = t_grid, y = x_dis, fdParobj = gamma_basis)

  # FPCA
  x_recv_pca  = fda::pca.fd(x_recv, nharm = num_pca)
  x_pca_score = x_recv_pca$scores

  # scale the x_pca_score.
  if(scale_kmean){
    scale_temp = scale(x_pca_score, center = FALSE)
    scale_vec  = attr(scale_temp, 'scaled:scale')
    x_s = scale_temp
  } else {
    x_s = x_pca_score
    scale_vec = rep(1, ncol(x_pca_score))
  }

  # If we center the observations subject-wise, then we do not have to estimate the ground mena mu and the fixed effects
  if(FALSE){
    temp_dat = data.frame(y, z, x_s, index[,3], ind = index[,1]) # index[,1] has the subject information.
    temp_dat = as.data.frame(temp_dat %>% group_by(ind) %>% mutate_all(.funs = scale, scale = FALSE))
    y_c = temp_dat[,1]
    z_c = as.matrix(temp_dat[,2:(1+ncol(z))])
    x_s_c = as.matrix(temp_dat[, (2+ncol(z)):(1+ ncol(z) + ncol(x_s))])
    treat_ind_c = temp_dat[, (2+ ncol(z) + ncol(x_s))]
  } else {
    y_c = y
    z_c = z
    x_s_c = x_s
    treat_ind_c = index[,3]
  }

  # a subject-wise summary matrix
  prodKron   = matrix(0,nrow=nrow(index),ncol=n_b)
  index_iter = 0
  struT      = (index %>% dplyr::group_by(ind_b) %>% dplyr::summarise(n_w=n()))$n_w
  for(i in 1:n_b ){
    prodKron[(index_iter+1):(index_iter+struT[i]),i] = 1;
    index_iter = index_iter + struT[i];
  }

  prodKron_mean   = prodKron %*% diag(1/struT)


  if(est_treat){
    z_treat_c = cbind(treat_ind_c, z_c)
    z_treat   = cbind(index[,3], z)
  } else {
    z_treat_c = z_c
    z_treat   = z
  }

  des_dummy       = create_design(index, est_fix_eff=est_fix_eff, mu = FALSE, est_treat = FALSE) # pseudo design, include fixed effects and treatment effects
  if(is.null(des_dummy)){
    des_dummy = create_design(index, est_fix_eff=TRUE, mu = FALSE, est_treat = FALSE)
    dummyhat  = rep(0, n_b)
  } else {
    dummyhat  = rep(0, ncol(des_dummy))  # the fixed effects only
  }

  fixeffhat_best = est_res$fixeffhat
  treathat_best  = est_res$treathat
  zhat_best      = est_res$zhat
  group_index_best  = est_res$group
  xhat_subject_best = est_res$xhat_subject
  obj_val_best  = est_res$aa # devided by sample size


  # systemetically relocate subjects to decrease object function
  uni_group = 1:num_group
  searching = TRUE
  num_search = 1
  while(searching & num_search < max_iter){
    searching = FALSE
    for(iter in 1:n_b){
      for(iter_g in setdiff(uni_group, group_index_best[iter])){
        # if move subject iter to another group, calculate the change in objective value.
        group_index_temp = group_index_best
        group_index_temp[iter] = iter_g

        x_score_group   = x_pca_m(x_s_c, group_index_temp, index) # for functional covariates
        X               = cbind(des_dummy, z_treat_c, x_score_group)    # pseduo design matrix

        names_dummy     = paste('dummy', 1:ncol(des_dummy), sep = '')
        names_z         = paste('z', 1:ncol(z_treat_c), sep = '')
        names_fun       = paste('fun', 1:ncol(x_score_group), sep = '')

        colnames(X)     = c(names_dummy, names_z, names_fun)

        lmfit           = lm(y_c ~ X + 0)       # the intercept is added automatically.
        coef_fit        = coef(lmfit)
        coef_fit[which(is.na(coef_fit))] = 0
        z_treat_hat     = coef_fit[paste('X',names_z, sep = '')]
        #xhat_group_old  = xhat_group
        xhat_group      = diag(1/scale_vec) %*% matrix(coef_fit[paste('X',names_fun, sep = '')], nrow=num_pca)
        if(est_fix_eff){
          #dummyhat        = c(t(prodKron_mean) %*% as.vector(y) - (t(prodKron_mean) %*% x_pca_m(x_pca_score, group_index_temp, index) ) %*% c(xhat_group) - (t(prodKron_mean) %*% z_treat) %*% z_treat_hat)   # the fixed effects estimate
          dummyhat = coef_fit[paste('X',names_dummy, sep = '')]
        } else {
          dummyhat = rep(0, n_b)
        }

        # the value of the object function
        #obj_val = sum((as.vector(y) - as.vector(des_dummy %*% dummyhat) - as.numeric(as.matrix(z_treat %*% z_treat_hat)) - x_pca_m(x_pca_score, group_index_temp, index) %*% c(xhat_group))^2) / length(y)
        obj_val  = log(sum(lmfit$residuals^2))
        if(obj_val < obj_val_best){
          uni_group_index  = (unique(group_index_temp))
          #xhat_subject_old = xhat_subject
          xhat_subject = xhat_subject_best
          for(ix in 1:n_b){
            xhat_subject[,ix] = xhat_group[, which(group_index_temp[ix] == uni_group_index)]
          }
          fixeffhat_best = dummyhat
          if(est_treat){
            treathat_best  = z_treat_hat[1]
            zhat_best      = z_treat_hat[2:length(z_treat_hat)]
          } else {
            treathat_best  = 0
            zhat_best      = z_treat_hat
          }
          group_index_best  = group_index_temp
          xhat_subject_best = xhat_subject
          obj_val_best  = obj_val # devided by sample size
          searching = TRUE
          break
        }
      }
    }
    num_search = num_search + 1
  }
  xhat_sub_fd   = fd(coef=x_recv_pca$harmonics$coefs %*% xhat_subject_best, basisobj = x_recv_pca$harmonics$basis)
  return(list(muhat = est_res$muhat, fixeffhat = fixeffhat_best, treathat = treathat_best, zhat=zhat_best,
              group=group_index_best, xhat_subject = xhat_subject_best, xhat_sub_fd=xhat_sub_fd, aa=obj_val_best,num_search=num_search))




}

#' Inference method for the partially functional linear model.
#' 
#' This function outputs the standard error of the scarlar estimates and the confidence band of the functional coefficient.
#' 
#' @param data_list A list of data. Several elements must be present in the list. The reponse \code{y},
#' the functional covariate \code{x_recv}, the scalar covariates \code{z}, and an index matrix \code{index}. 
#' The functional covariate \code{x_recv} must be generated from the \code{fda} package by, e.g., spline smoothing. 
#' The scalar covariates \code{z} is a matrix. The index matrix \code{index} is a data.frame recording the structure of the
#' data. The first column of \code{index} is the family number, the second column is the within family index. The column names of \code{index}
#' must be \code{ind_b} and \code{ind_w}.
#' 
#' @param group The estimated group membership from the function \code{bic_kmean_est}.
#' 
#' @param num_pca The number of principal components.
#' 
#' @param tau_1 The confidence level is \code{1-tau_1}, can be a vector.
#' 
#' @param tau_2 The proportion of the index of the functional coefficient that might be not covered by the band.
#' 
#' @param boot_size The size of bootstrap sampling.
#' 
#' 
#' 
#' @export
kmean_infer = function(data_list, group, num_pca = 5,  tau_1 = c(0.005, 0.05, 0.1), tau_2 = 0.1, boot_size = 5000){
  # inference for the parameters given a group structure.
  # One needs certain amount of undersmooth to reduce the bias in when construct confidence intervals/bands.
  y      = data_list$y
  #x_dis  = data_list$x_dis
  x_recv = data_list$x_recv # a functional object
  #t_grid = data_list$t_grid
  z      = data_list$z
  index  = data_list$index
  n_b    = length(unique(index[,1]))
  n      = nrow(z)

  est_treat = FALSE
  
  # FPCA
  x_recv_pca  = fda::pca.fd(x_recv, nharm = num_pca)
  x_pca_score = x_recv_pca$scores

  # center the data subject-wise to remove fixed effects.
  subj_center = TRUE
  if(subj_center){
    temp_dat = data.frame(y, z, x_pca_score, index[,3], ind = index[,1]) # index[,1] has the subject information.
    temp_dat = as.data.frame(temp_dat %>% group_by(ind) %>% mutate_all(.funs = scale, scale = FALSE))
    y = temp_dat[,1]
    z = as.matrix(temp_dat[,2:(1+ncol(z))])
    x_pca_score = as.matrix(temp_dat[, (2+ncol(z)):(1+ ncol(z)+ ncol(x_pca_score))])
    index[,3]   = temp_dat[, (2+ ncol(z) + ncol(x_pca_score))]
  }

  des_dummy     = create_design(index, est_fix_eff=FALSE, mu = FALSE, est_treat = est_treat) # pseudo design, include fixed effects and treatment effects

  # estimation
  x_score_group   = x_pca_m(x_pca_score, group, index) # for functional covariates
  X               = cbind(des_dummy, z, x_score_group)    # pseduo design matrix
  lmfit           = lm(y ~ X + 0)       # the intercept is added automatically.
  coef_fit        = coef(lmfit)
  coef_fit[which(is.na(coef_fit))] = 0

  # residual
  res_fit = residuals(lmfit)

  # bootstrap samples
  err_boot_m = matrix(sample(res_fit, n*boot_size, replace = TRUE), nrow = n, ncol = boot_size)

  # QR decomposition and hat matrix.
  qr_f  = base::qr(x_score_group)
  hat_f = base::tcrossprod(qr.Q(qr_f))

  # there are two parts
  z_treat = cbind(des_dummy, z)
  rhs     = err_boot_m - z_treat %*% solve(t(z_treat) %*% (diag(n) - hat_f) %*% z_treat, t(z_treat) %*% (diag(n) - hat_f) %*% err_boot_m )
  len_i   = x_recv$basis$rangeval[2] - x_recv$basis$rangeval[1]
  a_vec   = qr.solve(qr_f, rhs) # should be of size m*boot_size

  # first scheme, mimic Kato (2018) but include covariates
  a_vec2  = a_vec^2
  # calculate bound for each group
  bnd_fun = c()

  for(ii in 1: length(unique(group))){
    bnd_fun = rbind(bnd_fun, sqrt( quantile(apply(a_vec2[((ii-1)*num_pca+1):( ii*num_pca ),], 2, sum), probs = 1- tau_1) / ( tau_2 * len_i) ))
  }

  # the asymptotic variance matrix of the finite dimensional coefficients. It is the inverse of the following matrix
  Bhat = t(z_treat) %*% (diag(n) - hat_f) %*% z_treat

  # or by bootstrap again to get the variance.
  rhs_gamma = solve( Bhat, t(z_treat) %*% (diag(n) - hat_f) %*% err_boot_m)
  cov_gamma_boot = cov(t(rhs_gamma)) # take square root and multiply by 1.96 to get a bound

  # or use quantiles directly
  quan_gamma = apply(rhs_gamma, 1, quantile, probs = c(0.005, 0.025, 0.05, 0.995, 0.975, 0.95))

  # second scheme, pointwise interval
  # t_grid = seq(x_recv$basis$rangeval[1], x_recv$basis$rangeval[2], len=n_grid)
  # efun_val = eval.fd(t_grid, x_recv_pca$harmonics)
  #
  # bnd_pw = c()
  # for(ii in 1: length(unique(group))){
  #   bnd_pw = rbind(bnd_pw, apply( efun_val %*% a_vec[((ii-1)*num_pca+1):( ii*num_pca ),], 1, quantile , probs = c(0.005, 0.025, 0.05, 0.995, 0.975, 0.95) ) )
  # }

  return(list(bnd_fun = bnd_fun, Bhat = Bhat, cov_gamma_boot = cov_gamma_boot, quan_gamma = quan_gamma))
}


fun_dist <- function(fun1, fun2){
  n_fun1 = ncol(fun1$coefs)
  n_fun2 = ncol(fun2$coefs)
  dist_m = matrix(0, nrow = n_fun1, ncol = n_fun2)
  for(i in 1:n_fun1){
    for(j in 1:n_fun2){
      dist_m[i, j] = fda::inprod(fun1[i] - fun2[j],  fun1[i] - fun2[j])
    }
  }
  return(dist_m)
}








