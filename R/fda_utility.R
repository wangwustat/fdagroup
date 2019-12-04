

my_randindex <- function(IG, resGroup){
  N <- length(IG)
  trueP <- 0; trueN <- 0;
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      if(IG[i]==IG[j] && resGroup[i]==resGroup[j]){
        trueP=trueP + 1;}
      if(IG[i]!=IG[j] && resGroup[i]!=resGroup[j]){
        trueN=trueN+1;}
    }
  }
  randIndex=2*(trueP+trueN)/(N*(N-1));
  return(randIndex)
}

#' @export
data_gen_fdagroup <- function(index, nbasis = 49, n_grid=100, dev = 0.5, sd_err=1,
                              mu = 0, fix_eff=NULL, treat_eff = -1, beta = c(2, -1, 1.5, 5, -1.7),
                              group_member=NULL, coef_g1 = NULL, coef_g2 = NULL, sigma_seq = NULL){
  # index is a two column matrix, indicating the index system of longitudinal data.
  # first colunn indexed between, second within, third column indicator of treatment. 1-treated, 0-untreated.
  if(dim(index)[2] != 3){
    stop('index should has three columns')
  }
  index_df  = as.data.frame(index)
  colnames(index_df) = c('ind_b', 'ind_w', 'ind_t')
  index_df  = dplyr::arrange(index_df, ind_b, ind_w)
  index[,1] = dplyr::dense_rank(index[,1]) # such that index of between subject are 1,2,...n_b
  n         = nrow(index_df)                        # total number of obs
  n_b       = length(unique(index_df[,1]))          # number of subjects

  # generating functional covariates. only use cosine functions. and mean is zero.
  gamma_basis = fda::create.fourier.basis(nbasis = 2*nbasis+1, dropind = c(1, seq(2, 2*nbasis, by=2)))
  sigma_seq_t = 1 / ( ((1:nbasis) - 0.5)* pi )
  if(is.null(sigma_seq)){
    sigma_seq   = sigma_seq_t # root of the eigenvalue. since it is used in rnorm(sd=sigma_seq)
  } else {
    sigma_seq = c(sigma_seq, 2 * sigma_seq_t[-c(1:length(sigma_seq))], rep(1e-14, 10))[1:nbasis]
  }
  x = fChange::fun_IID(n, nbasis, Sigma = sigma_seq, basis = gamma_basis)

  # discretize x and reconstruct using fourier basis.
  t_grid  = seq(0, 1, len=n_grid)   # sampling grid
  x_dis   = fda::eval.fd(t_grid, x)  # sampling value
  x_basis = fda::create.fourier.basis(nbasis = 2*40 + 1, dropind = c(1, seq(2, 2*40, by=2)))
  #x_recv  = fda::smooth.basis(argvals = t_grid, y = x_dis, fdParobj = x_basis)$fd
  x_recv  = x

  # create the true functional coefficient. only 2 groups.
  # fda:fd(coef=NULL) coef can be a matrix,

  if(is.null(coef_g1)){
    #cvecf_1     = matrix(2^1.5*((-1)^(1:nbasis))*((1:nbasis)^(-2)), nbasis, 1)
    cvecf_1     =matrix(2^1.5*((-1)^(1:nbasis))*((1:nbasis+1)^(-2)), nbasis, 1) * 3
  } else {
    cvecf_1     = matrix(c(coef_g1, rep(0, nbasis))[1:nbasis], nbasis, 1)
  }
  #gamma_fd_1    = fda::fd(cvecf_1, gamma_basis)

  if(is.null(coef_g2)){
    #cvecf_2_t   = matrix(2*((-1)^(1:nbasis))*(sqrt(3)*(1:nbasis)^(-4) - sqrt(5)*(1:nbasis)^(-2)), nbasis, 1)
    cvecf_2_t   =   matrix(5^1.5*((-1)^(1:nbasis+1))*((1:nbasis+1)^(-3)), nbasis, 1) * 3
  } else {
    cvecf_2_t   = matrix(c(coef_g2, rep(0, nbasis))[1:nbasis], nbasis, 1)
  }
  cvecf_2     = cvecf_1 %*% (1-dev) + cvecf_2_t %*% dev
  #gamma_fd_2  = fda::fd(cvecf_2, gamma_basis)
  gamma_fd    = fda::fd(cbind(cvecf_1, cvecf_2), gamma_basis)

  # some more covariates and coefficients. z should be correlated with x, or
  # you can do regression analysis separately.
  p_z      = length(beta)
  zx_basis = fda::create.fourier.basis(nbasis = 11)
  zx_coeff = abs(matrix(rnorm(zx_basis$nbasis * p_z, sd = 1), zx_basis$nbasis, p_z))
  zx_fd    = fda::fd(zx_coeff, zx_basis)
  #z        = matrix(rnorm(n*p_z), n, p_z) * abs(fda::inprod(x, zx_fd))
  z        = matrix(rnorm(n*p_z), n, p_z)
  #beta     = c(2, -1, 1.5, 5, -1.7)

  # generate the response.
  error       = rnorm(n, sd = sd_err)
  if( is.null(group_member) ){
    # number of groups = 1+length(dev)
    group_member = adjust_group(sample(1:(1+length(dev)), n_b, replace = TRUE))
  }
  if(is.null(fix_eff)){
    fix_eff   = rnorm(n_b, sd = 1)
    fix_eff[1]= 0
  }
  #
  rand_eff_b  = data.frame(unique(index_df[,1]), fix_eff, as.factor(group_member)) # for each subject, generate a random effect and the group index
  colnames(rand_eff_b) = c('ind_b', 'fe', 'group')
  fe_eff      = dplyr::left_join(index_df, rand_eff_b, by='ind_b')$fe # fixed effects
  group_w     = dplyr::left_join(index_df, rand_eff_b, by='ind_b')$group
  group_w_m   = model.matrix(~group_w+0)

  #y           = mu + fe_eff + treat_eff*index_df[,3] + z %*% beta +  two_group * fda::inprod(x, gamma_fd_1) + (1-two_group) * fda::inprod(x, gamma_fd_2) + error
  y           = mu + fe_eff + treat_eff*index_df[,3] + z %*% beta +  apply(fda::inprod(x, gamma_fd) * group_w_m, 1, sum) + error
  return(list(y=y, x=x, z=z, x_recv=x_recv, t_grid=t_grid, index = index_df, group=group_member,
              gamma_fd=gamma_fd, mu=mu, beta=beta, fix_eff=fix_eff, treat_eff=treat_eff))
}

data_processing <- function(data_list){
  # adjust the index, such that it has a standardized form.
  # the between subject index is 1, 2, ..., n_b
  # within each subject, the within subject index is 1, 2, ..., n_w
  # the value of the treatment effect is 1-treated, 0-placebo. Only these two values are allowed currently.
  # finally, index should be a data.frame.
  index  = data_list$index

  # if index only has two columns, add one columns of 1s
  if(ncol(index) ==1 ){
    stop('index has too few columns (1)')
  } else if(ncol(index) == 2){
    index[, 3] = 1
  }

  colnames(index)[1:3] = c('ind_b', 'ind_wt', 'ind_t')
  index  = dplyr::arrange(index, ind_b, ind_wt)
  index[,1] = dplyr::dense_rank(index[,1]) # such that index of between subject are 1,2,...n_b
  index  = index %>% dplyr::group_by(ind_b) %>% dplyr::mutate(ind_w=1:length(ind_b)) %>% dplyr::select(ind_b, ind_w, ind_t)

  # if(length(unique(index$ind_t)) != 2){
  #   stop("more than two treatment groups are not allowed.")
  # }
  data_list$index = as.data.frame(index)
  return(data_list)
}


create_design = function(index, est_fix_eff = TRUE, mu = TRUE, est_treat = TRUE){
  res   = c()
  n_b   = length(unique(index[,1]))
  if(mu){
    res = cbind(res, matrix(1, nrow = nrow(index), ncol=1))
  }
  # first, fixed effect
  if(est_fix_eff){
    prodKron   = matrix(0,nrow=nrow(index),ncol=n_b)
    index_iter = 0
    struT      = (index %>% dplyr::group_by(ind_b) %>% dplyr::summarise(n_w=n()))$n_w
    for(i in 1:n_b ){
      prodKron[(index_iter+1):(index_iter+struT[i]),i] = 1;
      index_iter = index_iter + struT[i];
    }
    if(mu){
      res = cbind(res, prodKron[,-1]) # delete the first column such that mean is used as the base to compare
    } else {
      res = cbind(res, prodKron)
    }
  }

  # second, treatment effect
  if(est_treat){
    res = cbind(res, index$ind_t)
  }

  return(res)
}

#' @export
adjust_group <- function(group_index){
  # the group index might be like this: 1 1 1 1 3 3 1 2 2 2 1 1 1 3 2 1 2 1 1 2
  # adjust the index such that the first group is represented by 1, the second group is represented by 2.
  temp_group = as.numeric(rep(1, length(group_index)))
  uni_group  = unique(group_index)
  for (i in 1:length(uni_group) ){
    temp_group[group_index == uni_group[i]] = i
  }
  return(as.numeric(temp_group))
}


create_fixeff_design <- function(struT){
  # This function creates a pseudo design matrix for random effects.
  # struT contains the number of observations for each suject.
  prodKron   = matrix(0, nrow=sum(struT), ncol=length(struT))
  index_iter = 0
  for(i in 1:length(struT) ){
    prodKron[(index_iter+1):(index_iter+struT[i]),i] = 1;
    index_iter = index_iter + struT[i];
  }
  return(prodKron)
}

# estimation


# folds <- split(sample(nrow(data), nrow(data),replace=FALSE), as.factor(1:K))


kfold_cv <- function(data_list, K=5, num_group = 2, num_pca = 5, nbasis = 20, gamma_basis=NULL){
  # Do k-fold cv, and return the rmse.
  y      = data_list$y
  #x_dis  = data_list$x_dis
  x_recv = data_list$x_recv
  #t_grid = data_list$t_grid
  z      = data_list$z
  index  = data_list$index

  n_b    = length(unique(index[,1]))
  struT  = (index %>% dplyr::group_by(ind_b) %>% dplyr::summarise(n_w=n()))$n_w

  # first generate partitions. The partition is done within subject. stratified sampling.
  folds  = vector('list', n_b)
  for(i in 1:n_b){
    folds_temp = split(sample(struT[i], struT[i],replace=FALSE), as.factor(1:K))
    folds[[i]] = folds_temp
  }

  # estimation.
  pred_error = c()
  for(iter_k in 1:K){ # iterate in K-folds
    # extract index
    index_est = c()
    index_pre = c()
    # extract sample index
    for(i_b in 1:n_b){
      temp = folds[[i_b]][[iter_k]]
      index_pre = rbind(index_pre,cbind(rep(i_b, length(temp)), temp))
    }
    index_pre = data.frame(index_pre)
    colnames(index_pre) = c('ind_b', 'ind_w')
    index_pre = dplyr::arrange(index_pre, ind_b, ind_w)
    index_est = dplyr::anti_join(index, index_pre, by = c("ind_b", "ind_w"))

    #extract data.
    data_kcv      = extract_data(data_list, index_pre, index_est)
    data_list_pre = data_kcv$pre
    data_list_est = data_kcv$est

    # estimation
    kcv_est = fdagroup_est(data_list_est, num_group = num_group, num_pca = num_pca, nbasis = nbasis, gamma_basis=gamma_basis)

    # prediction
    # literally first column of xhat_group is the coefficient for the group '1', second for group '2', etc,...
    pred_obj = predict_fdagroup(kcv_est, data_list_pre)
    pred_error = c(pred_error, pred_obj$pred_error)
  }
  # K-fold CV average prediction error.
  return(mean(pred_error^2))
}


kfold_compare <- function(data_list, K=5, num_group = 2, num_pca = 5, est = NULL){
  # Do k-fold cv, and return the rmse.
  # May be we can treat the estimated groups as fixed.
  # compare the K-mean model, the pooled model, and the family-wise model on the ability of prediction.
  y      = data_list$y
  #x_dis  = data_list$x_dis
  x_recv = data_list$x_recv
  #t_grid = data_list$t_grid
  z      = data_list$z
  index  = data_list$index

  n_b    = length(unique(index[,1]))
  struT  = (index %>% dplyr::group_by(ind_b) %>% dplyr::summarise(n_w=n()))$n_w

  # first generate partitions. The partition is done within subject. stratified sampling.
  folds  = vector('list', n_b)
  for(i in 1:n_b){
    folds_temp = split(sample(struT[i], struT[i],replace=FALSE), as.factor(1:K))
    folds[[i]] = folds_temp
  }

  # estimation.
  pred_error = c()
  for(iter_k in 1:K){ # iterate in K-folds
    # extract index
    index_est = c()
    index_pre = c()
    # extract sample index
    for(i_b in 1:n_b){
      temp = folds[[i_b]][[iter_k]]
      index_pre = rbind(index_pre,cbind(rep(i_b, length(temp)), temp))
    }
    index_pre = data.frame(index_pre)
    colnames(index_pre)[1:2] = c('ind_b', 'ind_w')
    index_pre = dplyr::arrange(index_pre, ind_b, ind_w)
    index_est = dplyr::anti_join(index, index_pre, by = c("ind_b", "ind_w"))

    #extract data.
    data_kcv      = extract_data(data_list, index_pre, index_est)
    data_list_pre = data_kcv$pre
    data_list_est = data_kcv$est

    # estimation
    if(is.null(est)){
      est_kmean_cv = bic_kmean_est(data_list_est, num_group = num_group, num_pca = num_pca, est_fix_eff=TRUE, est_treat = FALSE)
    } else { # a pre-estimate has been given, fix the estimated groups.
      est_kmean_cv = bic_kmean_est(data_list_est, num_group = num_group, num_pca = num_pca, est_fix_eff=TRUE, est_treat = FALSE, group_index = est$group, max_iter = 1)
    }
    est_subje_cv   = sel_subjectwise(data_list_est, num_pca_vec = num_pca, est_fix_eff = TRUE, est_treat = FALSE)
    est_pool_cv    = bic_kmean_est(data_list_est, num_group = 1, num_pca = num_pca, est_treat = FALSE)

    # prediction
    fit_obj_kmean  = predict_fdagroup(est_kmean_cv$est8, data_list_pre)
    fit_obj_family = predict_fdagroup(est_subje_cv$est8, data_list_pre)
    fit_obj_pool   = predict_fdagroup(est_pool_cv$est8, data_list_pre)

    # literally first column of xhat_group is the coefficient for the group '1', second for group '2', etc,...
    perror_temp = cbind(index_pre, fit_obj_kmean$pred_error, fit_obj_family$pred_error, fit_obj_pool$pred_error)
    colnames(perror_temp) = c('ind_b', 'ind_w', 'perr_kmean', 'perr_family', 'perr_pool')
    pred_error = rbind(pred_error, perror_temp)
  }
  pred_error = dplyr::arrange(pred_error, ind_b, ind_w)
  # K-fold CV average prediction error.

  return(pred_error)
}


cv_compare <- function(data_list, K=5, num_group = 2, num_pca = 5, est = NULL){
  # Do a traditional leave one out cv
  # compare the K-mean model, the pooled model, and the family-wise model on the ability of prediction.
  y      = data_list$y
  #x_dis  = data_list$x_dis
  x_recv = data_list$x_recv
  #t_grid = data_list$t_grid
  z      = data_list$z
  index  = data_list$index

  n_b    = length(unique(index[,1]))
  struT  = (index %>% dplyr::group_by(ind_b) %>% dplyr::summarise(n_w=n()))$n_w

  # estimation.
  pred_error = c()
  for(iter_k in 1:nrow(index)){ # iterate in K-folds
    # extract sample index
    index_pre = index[iter_k,,drop=FALSE]
    index_est = index[-iter_k,]

    #extract data.
    data_kcv      = extract_data(data_list, index_pre, index_est)
    data_list_pre = data_kcv$pre
    data_list_est = data_kcv$est

    # estimation
    if(is.null(est)){
      est_kmean_cv = bic_kmean_est(data_list_est, num_group = num_group, num_pca = num_pca, est_fix_eff=TRUE, est_treat = FALSE)
    } else {
      est_kmean_cv = bic_kmean_est(data_list_est, num_group = num_group, num_pca = num_pca, est_fix_eff=TRUE, est_treat = FALSE, group_index = est$group, max_iter = 1)
    }
    est_subje_cv = sel_subjectwise(data_list_est, num_pca_vec = num_pca, est_fix_eff = TRUE, est_treat = FALSE)
    est_pool_cv  = bic_kmean_est(data_list_est, num_group = 1, num_pca = num_pca, est_treat = FALSE)

    # prediction
    fit_obj_kmean  = predict_fdagroup(est_kmean_cv$est8, data_list_pre)
    fit_obj_family = predict_fdagroup(est_subje_cv$est8, data_list_pre)
    fit_obj_pool   = predict_fdagroup(est_pool_cv$est8, data_list_pre)

    # literally first column of xhat_group is the coefficient for the group '1', second for group '2', etc,...
    perror_temp = cbind(index_pre, fit_obj_kmean$pred_error, fit_obj_family$pred_error, fit_obj_pool$pred_error)
    colnames(perror_temp) = c('ind_b', 'ind_w', 'ind_t', 'perr_kmean', 'perr_family', 'perr_pool')
    #perror_temp = cbind(index_pre, fit_obj_kmean$pred_error, fit_obj_pool$pred_error)
    #colnames(perror_temp) = c('ind_b', 'ind_w', 'perr_kmean', 'perr_pool')
    pred_error = rbind(pred_error, perror_temp)
  }
  #pred_error = dplyr::arrange(pred_error, ind_b, ind_w)
  # K-fold CV average prediction error.

  return(pred_error)
}

# crossvalidation and estimation.

cv_est <- function(data_list, K=5, num_group = 2, num_pca = 5, nbasis = 20, gamma_basis=NULL){
  # num_group and num_pca could be a grid, by k-cv, select the best combination and return estmated coefficients.
  #cat("In the function cv_est\n")
  cv_error = c()
  para_df  = expand.grid(num_group, num_pca)
  #saveRDS(data_list, 'err_data.rds')
  for( iter in 1:nrow(para_df)){
    num_group_t = para_df[iter, 1]
    num_pca_t   = para_df[iter, 2]
    #print(sprintf('num_group %d num_pca %d', num_group_t, num_pca_t))
    cv_temp     = kfold_cv(data_list, K=K, num_group = num_group_t, num_pca = num_pca_t, nbasis = nbasis, gamma_basis=gamma_basis)
    cv_error    = c(cv_error, cv_temp)
  }

  # the parameter combination which minimize the cv-error.
  cv_index    = which.min(cv_error)
  num_group_t = para_df[cv_index, 1]
  num_pca_t   = para_df[cv_index, 2]
  cv_est      = fdagroup_est(data_list, num_group = num_group_t, num_pca = num_pca_t, nbasis = nbasis, gamma_basis=gamma_basis)
  return(list(cv_est=cv_est, num_group=num_group_t, num_pca=num_pca_t))
}

# prediction

predict_fdagroup <- function(kcv_est, data_pre, err = TRUE){
  if(err){ # whether to compute the prediction error
    y    = data_pre$y
  }
  #x_dis  = data_pre$x_dis
  x_recv = data_pre$x_recv
  mn_fd  = mean(x_recv)
  #x_recv = x_recv - mean(x_recv) # centralize, because it is centered when doing PCA
  #t_grid = data_pre$t_grid
  z      = data_pre$z
  #gamma_basis = kcv_est$gamma_basis
  index  = data_pre$index
  colnames(index) = c('ind_b', 'ind_w', 'ind_t')
  #index_df  = dplyr::arrange(index, ind_b, ind_w)
  index_df  = data.frame(index, tt_ind=1:nrow(index))
  sub_uni   = unique(index[,1])
  struT     = (index %>% dplyr::group_by(ind_b) %>% dplyr::summarise(n_w=n()))$n_w

  # literally first column of xhat_group is the coefficient for the group '1', second for group '2', etc,...
  group_hat  = kcv_est$group
  pred_error = c()
  yfit       = c()
  func_fit   = c()
  for(i_b in sub_uni){
    # select sample in the i_b-th subject
    ind_wsub = dplyr::filter(index_df, ind_b == i_b)$tt_ind
    if(length(ind_wsub) == 0){
      next
    }
    if(err){
      y_wsub = y[ind_wsub]
    }
    #x_dis_wsub = x_dis[, ind_wsub,drop=FALSE]
    z_wsub = z[ind_wsub,,drop=FALSE]
    x_recv_wsub = x_recv[ind_wsub]
    treat_wsub  = index[ind_wsub, 3, drop=FALSE]

    # smooth x
    #x_recv  = fda::smooth.basis(argvals = t_grid, y = x_dis_wsub, fdParobj = gamma_basis)

    # fit
    yfit_temp  = kcv_est$muhat + kcv_est$fixeffhat[i_b] + c(as.matrix(treat_wsub) %*% kcv_est$treathat) + as.numeric(z_wsub %*% kcv_est$zhat) +
           c(fda::inprod(x_recv_wsub, kcv_est$xhat_sub_fd[i_b])) - c(fda::inprod(mn_fd, kcv_est$xhat_sub_fd[i_b]))
    if(err){
      pred_error = c(pred_error, y_wsub-yfit_temp)
    }
    yfit       = c(yfit, yfit_temp)
    func_fit   = c(func_fit, c(fda::inprod(x_recv_wsub, kcv_est$xhat_sub_fd[i_b])) )
  }
  return(list(yfit = yfit, pred_error=pred_error, index=index, func_fit = func_fit))
}


# extract data in the K-fold cv

extract_data <- function(data_list, index_pre, index_est){
  y      = data_list$y
  #x_dis  = data_list$x_dis
  x_recv = data_list$x_recv
  #t_grid = data_list$t_grid
  z      = data_list$z
  index  = data_list$index
  #colnames(index) = c('ind_b', 'ind_w')
  index_df  = dplyr::arrange(index, ind_b, ind_w)
  index_df  = data.frame(index_df, tt_ind=1:nrow(index))

  tt_pre = dplyr::semi_join(index_df, index_pre, by=c('ind_b', 'ind_w'))$tt_ind
  tt_est = dplyr::semi_join(index_df, index_est, by=c('ind_b', 'ind_w'))$tt_ind

  y_pre  = y[tt_pre]
  #x_dis_pre = x_dis[, tt_pre] # columns of x are discretized observations of functions
  z_pre  = z[tt_pre,,drop = FALSE]
  x_recv_pre = x_recv[tt_pre]
  index_pre_f = index[tt_pre,,drop=FALSE]

  y_est  = y[tt_est]
  #x_dis_est = x_dis[, tt_est]
  z_est  = z[tt_est,,drop=FALSE]
  x_recv_est = x_recv[tt_est]
  index_est_f = index[tt_est,,drop=FALSE ]

  data_list_pre = list(y=y_pre, z=z_pre, x_recv=x_recv_pre, index = index_pre_f)
  data_list_est = list(y=y_est, z=z_est, x_recv=x_recv_est, index = index_est_f)

  return(list(est=data_list_est, pre=data_list_pre ))
}

# data generation and estimation
#' @export
datagen_est <- function(M=100, n_b=20, n_w=10, nbasis = 49, num_group_vec = 1:4,num_pca_vec = seq(2, 10, 5), lamGroup = c(1),
                        K=5, n_grid=100, dev = 0.5, sd_err=1, mu = 0, fix_eff=NULL, treat_eff = -1,
                        beta = c(2, -1, 1.5, 5, -1.7), group_member = group_member,
                        tuPara = c(0.1, 1, 3, 1, 10000, 2), precision = c(1e-5, 1e-5),
                        coef_g1 = c(2, 1, 3, 4), coef_g2 = c(-1, 2, -3, 4), des = "",
                        lam_max_interval = c(0.1, 10), lam_len = 100, sigma_seq = NULL){
  #index, p_z=5, nbasis = 49, n_grid=100, dev = 0.5, sd_err=1, sd_re = 0.5
  # index -- the subject and within subject index. p_z -- the number of covariate z. n_grid, how many discretized obs each function has.
  # dev -- indicate the dissimilarity of the two groups. sd_err -- error standard deviation. sd_re -- random effects standard deviation.

  # generate index of data.
  index = cbind(expand.grid(1:n_w, 1:n_b)[, c(2,1)], rep(c(rep(1, floor(n_w/2)), rep(0, n_w-floor(n_w/2))), n_b))
  colnames(index) = c('ind_b', 'ind_w', 'ind_t')  # index of between, within subjects and index of treatment.

  res <- foreach::foreach( iter = 1:M,.packages=c('fdacluster', 'fda', 'Matrix', 'tidyverse', 'Rcpp') ) %dopar% {
    #for(i in 1:M){
    data_list = data_gen_fdagroup(index, nbasis = nbasis, n_grid=n_grid, dev = dev, sd_err=sd_err,
                                  mu = mu, fix_eff=fix_eff, treat_eff = treat_eff, beta = beta, group_member = group_member,
                                  coef_g1 = coef_g1, coef_g2 = coef_g2, sigma_seq = sigma_seq)
    data_list = data_processing(data_list)
    my_time   = c()
    # estimation
    ptm <- proc.time()
    est_kmean = bic_kmean_est(data_list, num_group = num_group_vec, num_pca = num_pca_vec, loc_search = FALSE, max_iter = 100)
    cat(sprintf("Kmean time"))
    temp_time = proc.time() - ptm
    my_time = c(my_time, temp_time[3])

    ptm <- proc.time()
    est_subje = sel_subjectwise(data_list, num_pca_vec, est_treat = FALSE)
    temp_time = proc.time() - ptm
    my_time = c(my_time, temp_time[3])

    ptm <- proc.time()
    est_pool  = bic_kmean_est(data_list, num_group = 1, num_pca = num_pca_vec, est_treat = FALSE)
    cat(sprintf("Subject-wise and Pooled estimate time"))
    print(proc.time() - ptm)
    temp_time = proc.time() - ptm
    my_time = c(my_time, temp_time[3])

    ptm <- proc.time()
    est_naive_kmean = naive_kmeans(data_list, num_group = num_group_vec, num_pca = num_pca_vec, est_treat = FALSE)
    cat(sprintf("Two step kmeans"))
    print(proc.time() - ptm)
    temp_time = proc.time() - ptm
    my_time = c(my_time, temp_time[3])

    interval_kmean = kmean_infer(data_list, group = est_kmean$est8$group, num_pca = 4, est_treat = FALSE, tau_1 = c(0.005, 0.05, 0.1), tau_2 = 0.1, boot_size = 5000, n_grid = 200)

    ptm <- proc.time()
    tuPara[4] = 2
    est_pena_scad  = pena_est_fda_scale(data_list, num_pca_vec, tuPara, est_fix_eff=TRUE, precision = precision, scale_pena = FALSE, est_treat = FALSE, lamGroup_a = lamGroup)
    cat(sprintf("SCAD time"))
    print(proc.time() - ptm)
    temp_time = proc.time() - ptm
    my_time = c(my_time, temp_time[3])

    ptm <- proc.time()
    tuPara[4] = 1
    est_pena_lasso = pena_est_fda_scale(data_list, num_pca_vec, tuPara, est_fix_eff=TRUE, precision = precision, scale_pena = FALSE, est_treat = FALSE, lamGroup_a = lamGroup)
    cat(sprintf("LASSO time"))
    print(proc.time() - ptm)
    temp_time = proc.time() - ptm
    my_time = c(my_time, temp_time[3])

    #x=data_list$x, z=data_list$z,
    data_store = list(y=data_list$y, group = data_list$group, mu=data_list$mu, beta=data_list$beta, fix_eff=data_list$fix_eff,
                      treat_eff=data_list$treat_eff, gamma_fd = data_list$gamma_fd, dev = dev)
    list(data = data_store, est_kmean = est_kmean, est_naive_kmean = est_naive_kmean, est_pena_scad = est_pena_scad, est_pena_lasso = est_pena_lasso,
         est_subje = est_subje, est_pool = est_pool, interval_kmean = interval_kmean, my_time = my_time)
  }
  if(length(dev) == 1){
    filename = sprintf('fda_n_nb%d_nw%d_dev%2.2f_sd_err%2.1f_%s_%s.rds', n_b, n_w, dev, sd_err, stringr::str_flatten(num_pca_vec), des)
  } else {
    filename = sprintf('fda_n_nb%d_nw%d_G%d_sd_err%2.1f_%s_%s.rds', n_b, n_w, length(dev)+1, sd_err, stringr::str_flatten(num_pca_vec), des)
  }
  saveRDS(res, file = filename)
}

#(index, p_z=5, nbasis = 49, n_grid=100, dev = 0.5, sd_err=1, sd_re = 0.5)


# data generation and estimation

datagen_time <- function(M=100, n_b=20, n_w=10, nbasis = 49, num_group_vec = 1:4,num_pca_vec = seq(2, 10, 5), lamGroup = c(1),
                        K=5, n_grid=100, dev = 0.5, sd_err=1, mu = 0, fix_eff=NULL, treat_eff = -1,
                        beta = c(2, -1, 1.5, 5, -1.7), group_member = group_member,
                        tuPara = c(0.1, 1, 3, 1, 10000, 2), precision = c(1e-5, 1e-5),
                        coef_g1 = c(2, 1, 3, 4), coef_g2 = c(-1, 2, -3, 4), des = "",
                        lam_max_interval = c(0.1, 10), lam_len = 100, sigma_seq = NULL){
  #index, p_z=5, nbasis = 49, n_grid=100, dev = 0.5, sd_err=1, sd_re = 0.5
  # index -- the subject and within subject index. p_z -- the number of covariate z. n_grid, how many discretized obs each function has.
  # dev -- indicate the dissimilarity of the two groups. sd_err -- error standard deviation. sd_re -- random effects standard deviation.

  # generate index of data.
  index = cbind(expand.grid(1:n_w, 1:n_b)[, c(2,1)], rep(c(rep(1, floor(n_w/2)), rep(0, n_w-floor(n_w/2))), n_b))
  colnames(index) = c('ind_b', 'ind_w', 'ind_t')  # index of between, within subjects and index of treatment.

  res <- foreach::foreach( iter = 1:M,.packages=c('fdacluster', 'fda', 'Matrix', 'tidyverse', 'Rcpp') ) %dopar% {
    #for(i in 1:M){
    data_list = data_gen_fdagroup(index, nbasis = nbasis, n_grid=n_grid, dev = dev, sd_err=sd_err,
                                  mu = mu, fix_eff=fix_eff, treat_eff = treat_eff, beta = beta, group_member = group_member,
                                  coef_g1 = coef_g1, coef_g2 = coef_g2, sigma_seq = sigma_seq)
    data_list = data_processing(data_list)
    my_time   = c()
    # estimation


    ptm <- proc.time()
    est_subje = sel_subjectwise(data_list, num_pca_vec, est_treat = FALSE)
    temp_time = proc.time() - ptm
    my_time = c(my_time, temp_time[3])

    ptm <- proc.time()
    est_pool  = bic_kmean_est(data_list, num_group = 1, num_pca = num_pca_vec, est_treat = FALSE)
    cat(sprintf("Subject-wise and Pooled estimate time"))
    print(proc.time() - ptm)
    temp_time = proc.time() - ptm
    my_time = c(my_time, temp_time[3])


    #x=data_list$x, z=data_list$z,
    data_store = list(group = data_list$group, mu=data_list$mu, beta=data_list$beta, fix_eff=data_list$fix_eff,
                      treat_eff=data_list$treat_eff, gamma_fd = data_list$gamma_fd, dev = dev)
    list(data = data_store, my_time = my_time)
  }
  if(length(dev) == 1){
    filename = sprintf('time_nb%d_nw%d_dev%2.2f_sd_err%2.1f_%s_%s.rds', n_b, n_w, dev, sd_err, stringr::str_flatten(num_pca_vec), des)
  } else {
    filename = sprintf('time_nb%d_nw%d_G%d_sd_err%2.1f_%s_%s.rds', n_b, n_w, length(dev)+1, sd_err, stringr::str_flatten(num_pca_vec), des)
  }
  saveRDS(res, file = filename)
}



