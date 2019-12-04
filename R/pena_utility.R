
subjectwise_fda <- function(data_list, num_pca=2, est_fix_eff=TRUE, est_treat=FALSE){
  # estimate under the assumption that all the subject has different functional coefficients
  # if sub TRUE, estimate subject-wise. if sub FALSE, estimate the pooled model.
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

  # dummy variables for the fixed effects and treatment effects.
  des_dummy = create_design(index, est_fix_eff=est_fix_eff, mu=TRUE, est_treat = est_treat)
  Zplus     = cbind(des_dummy, z)   # this is the design matrix for all coefficients without groups.

  # create a diagonal matrix where each block consists of the observations of x_pca_score for each subject.
  struT  = (index %>% dplyr::group_by(ind_b) %>% dplyr::summarise(n_w=n()))$n_w
  diagX  = matrix(0, nrow=sum(struT), ncol=n_b*ncol(x_pca_score))
  countr = 0
  countl = 0
  for (i in 1:n_b){
    diagX[(1+countr):(countr+struT[i]),(1+countl):(countl+ncol(x_pca_score)) ] = x_pca_score[(1+countr):(countr+struT[i]),]
    countr = countr + struT[i]
    countl = countl + ncol(x_pca_score)
  }

  # estimation using the lm module
  X     = cbind(Zplus, diagX)
  lmfit = lm(y ~ X+0)
  bhat  = as.numeric(coef(lmfit))
  bhat[which(is.na(bhat))] = 0

  # gound mean, fixed effects and treatment effects
  group_index = 1:n_b
  res_list    = pack_pena(n_b, ncol(z), x_recv_pca, bhat, group_index, est_fix_eff=est_fix_eff, est_treat = est_treat)

  return(res_list)
}


#' Estimate a partially functional linear regression model with family-specific functional coefficients.
#' 
#' The estimator assumes family-specific functional coefficients.  
#' The the number of principal components are determined by the BIC criterion.
#' 
#' @param data_list A list of data. Several elements must be present in the list. The reponse \code{y},
#' the functional covariate \code{x_recv}, the scalar covariates \code{z}, and an index matrix \code{index}. 
#' The functional covariate \code{x_recv} must be generated from the \code{fda} package by, e.g., spline smoothing. 
#' The scalar covariates \code{z} is a matrix. The index matrix \code{index} is a data.frame recording the structure of the
#' data. The first column of \code{index} is the family number, the second column is the within family index. The column names of \code{index}
#' must be \code{ind_b} and \code{ind_w}.
#' 
#' @param num_pca A vector of candidate number of principal components.
#' 
#' @param est_fix_eff A logical value. If \code{TRUE}, then the fixed effects are estimated. Otherwise,
#' the fixed effects are not estimated
#' 
#' 
#' @export
family_wise <- function(data_list, num_pca_vec, est_fix_eff=TRUE){
  # select the number of pca in the subjectwise estimation by BIC
  #cat("In the function sel_subjectwise\n")
  bicvalue = rep(NA, length(num_pca_vec))
  n     = nrow(data_list$index) # total number of observations
  n_b   = length(unique(data_list$index[,1]))
  if(est_fix_eff){
    n_f = n_b + 1 + ncol(data_list$z)
  } else {
    n_f = 2 + ncol(data_list$z)
  }

  est1 = NULL; bicValue1 = 1e20; lambdabic1 = 1

  a_vector    = c()
  sumq_vector = c()
  b_vector    = c()

  est_treat = FALSE
    bic_factor  = 1
  

  for(num_pca in num_pca_vec){
    # estimate step
    subfit = subjectwise_fda(data_list, num_pca, est_fix_eff=est_fix_eff, est_treat = FALSE)

    # bic value
    aa   = mean(predict_fdagroup(subfit, data_list)$pred_error^2, na.rm = TRUE) # RSS
    sumq = num_pca * n_b + ncol(data_list$z) + as.numeric(est_fix_eff) * n_b + as.numeric(est_treat) # number of parameters.

    #bic_factor  = seq(from = 0.30, 1.0, length.out = 8)
    tbicValue   = aa + bic_factor *log(n) * (sumq)/n  #BIC

    a_vector    = c(a_vector, aa)
    b_vector    = c(b_vector, log(n) * (sumq)/n)
    sumq_vector = c(sumq_vector, sumq)

    for(i in 1:length(bic_factor)){
      if(tbicValue[i] < eval(parse(text = paste("bicValue", i, sep = "")))){
        assign(paste("est", i, sep = ""), subfit)
        assign(paste("bicValue", i, sep = ""), tbicValue[i])
        assign(paste("lambdabic", i, sep = ""), num_pca)
      }
    }
  }
  return(list(est1 = est1, bicValue1 = bicValue1, lambdabic1 = lambdabic1))
}

create_penamatrix <- function(n_b, p_x){
  # create a sparse matrix Delta.
  spi = rep(0,n_b*(n_b-1))
  spj = rep(0,n_b*(n_b-1))
  sps = rep(0,n_b*(n_b-1))
  count = 1
  index = 1
  for (i in 1:(n_b-1)){
    for (j in (i+1):n_b){
      spi[index]=count; spj[index]=i; sps[index]=1; index=index+1;
      spi[index]=count; spj[index]=j; sps[index]=-1; index=index+1;
      count=count+1;
    }
  }
  Delta = sparseMatrix(i=spi,j=spj,x=sps,dims=c(n_b*(n_b-1)/2,n_b), giveCsparse=TRUE)
  A = Delta %x% diag(nrow=p_x)
  return(list(A=A, Delta = Delta))
}

initial_pena <- function(data_list, num_pca, A, est_fix_eff=TRUE, est_treat=TRUE){
  est = subjectwise_fda(data_list, num_pca=num_pca, est_fix_eff=est_fix_eff, est_treat = est_treat)
  if(est_fix_eff){
    if(est_treat){
      theta = as.numeric(c(est$muhat, est$fixeffhat[-1], est$treathat, est$zhat, as.numeric(est$xhat_group)))
    }else{
      theta = as.numeric(c(est$muhat, est$fixeffhat[-1], est$zhat, as.numeric(est$xhat_group)))
    }
  } else {
    if(est_treat){
      theta = as.numeric(c(est$muhat, est$treathat, est$zhat, as.numeric(est$xhat_group)))
    } else {
      theta = as.numeric(c(est$muhat, est$zhat, as.numeric(est$xhat_group)))
    }
  }
  eta  = as.numeric(A %*% theta) + rnorm(nrow(A), sd=0.5)
  lagU = rep(0, nrow(A))
  return(list(theta = theta, eta = eta, lagU = lagU))
}

pack_pena <- function(n_b, p_z, x_recv_pca, theta, group_index, est_fix_eff=TRUE, est_treat=TRUE){
  # pack the results into a standardized form.
  # theta is vector include all the parameter estimates.
  muhat = theta[1]
  if(est_fix_eff){
    if(est_treat){
      fixeffhat = c(0, theta[2:n_b]) # now the fixed effects has length n_b, and the first element is 0
      treathat  = theta[n_b+1]
      zhat      = theta[(n_b+2):(n_b+1+p_z)]
      xhat_temp = matrix(theta[(n_b+2+p_z):length(theta)], ncol = n_b)
      xhat_sub_fd = fd(coef=x_recv_pca$harmonics$coefs %*% xhat_temp, basisobj = x_recv_pca$harmonics$basis)
    } else {
      fixeffhat = c(0, theta[2:n_b]) # now the fixed effects has length n_b, and the first element is 0
      treathat  = 0
      zhat      = theta[(n_b+1):(n_b+p_z)]
      xhat_temp = matrix(theta[(n_b+1+p_z):length(theta)], ncol = n_b)
      xhat_sub_fd = fd(coef=x_recv_pca$harmonics$coefs %*% xhat_temp, basisobj = x_recv_pca$harmonics$basis)
    }
  } else {
    if(est_treat){
      fixeffhat = rep(0, n_b)
      treathat  = theta[2]
      zhat      = theta[3:(p_z+2)]
      xhat_temp = matrix(theta[(p_z+3):length(theta)], ncol = n_b)
      xhat_sub_fd = fd(coef=x_recv_pca$harmonics$coefs %*% xhat_temp, basisobj = x_recv_pca$harmonics$basis)
    } else {
      fixeffhat = rep(0, n_b)
      treathat  = 0
      zhat      = theta[2:(p_z+1)]
      xhat_temp = matrix(theta[(p_z+2):length(theta)], ncol = n_b)
      xhat_sub_fd = fd(coef=x_recv_pca$harmonics$coefs %*% xhat_temp, basisobj = x_recv_pca$harmonics$basis)
    }
  }
  return(list(muhat = muhat, fixeffhat = fixeffhat, treathat = treathat, zhat=zhat,
              xhat_group=xhat_temp, xhat_sub_fd=xhat_sub_fd, group=group_index))
}

mynorm <- function(x){
  # 2-norm of x
  return(sqrt(crossprod(x,x)))
}

group_ind <- function(n_b, p_x, eta, Delta, thresh = 1e-4){
  # return the estimated group index.
  tempgroup = 1:n_b
  count =  1
  for ( ii in 1:(length(eta)/p_x) ){
    if (mynorm(eta[((count-1)*p_x+1):(count*p_x)]) < thresh){
      groupindex = sort(which(abs(Delta[ii,])>0.5))
      tempgroup[groupindex[2]] = tempgroup[groupindex[1]];
    }
    count=count+1;
  }
  return(tempgroup)
}


plot_bic <- function(bic_factor, a_vector, b_vector, lamGroup){
  #bic_factor <- seq(from = 0.35, 0.7, length.out = 8)
  bic_frame  = matrix(b_vector, ncol = 1) %*% matrix(bic_factor, nrow = 1) + a_vector
  colnames(bic_frame) = paste('BIC', bic_factor)
  bic_frame  = data.frame(bic_frame, lambda = lamGroup)

  bic_frame %>% tidyr::gather(key = 'bic', value = 'value', 1:8) %>% ggplot() + geom_point(aes(x=lambda, y=value)) + facet_wrap(~bic)

}


find_max_lam <- function(interval, n_b, num_pca, Delta,  y, Z, x_pca, struT, A, theta, eta, lagU, tuPara, precision){
  # find the upper limit of the lambda values.
  # upper limit of lambda should make 1 group
  # interval is to be searched.
  fun_max <- function(lam,n_b, num_pca, y, Z, x_pca, struT, A, theta, eta, lagU, tuPara, precision){
    tuPara[1] = lam
    bicRes    = fdagrouping(y, Z, x_pca, struT, A, theta, eta, lagU, tuPara, precision)
    eta       = bicRes$eta
    tempgroup = adjust_group(group_ind(n_b, num_pca, eta, Delta, thresh = 1e-4))
    return( abs(length(unique(tempgroup)) - 1) + lam/n_b )
  }
  lam_op = optimize(fun_max, interval = interval, tol = 1e-2, n_b = n_b, num_pca = num_pca, y=y, Z=Z, x_pca=x_pca, struT=struT,
                    A=A, theta=theta, eta=eta, lagU=lagU, tuPara=tuPara, precision=precision )
  return(lam_op$minimum)
}



pena_est_fda_path <- function(data_list, lamGroup, num_pca_vec, tuPara, est_fix_eff=TRUE, est_treat = TRUE, precision = c(1e-5, 1e-5),
                              bic_factor = NULL){
  y      = data_list$y
  #x_dis  = data_list$x_dis
  x_recv = data_list$x_recv
  #t_grid = data_list$t_grid
  z      = data_list$z
  index  = data_list$index
  n_b    = length(unique(index[,1]))
  struT  = as.numeric((index %>% dplyr::group_by(ind_b) %>% dplyr::summarise(n_w=n()))$n_w)
  n      = nrow(data_list$index) # total number of observations

  #cat("In the function pena_est_fda\n")
  #saveRDS(data_list, file='fail.rds')

  # create the pseudo design matrix
  des_dummy   = create_design(index, est_fix_eff=est_fix_eff, mu=TRUE, est_treat=est_treat) # pseudo design, include mean, fixed effects and treatment effects
  Z           = cbind(des_dummy, z)

  # store the estimators
  est1 = c();    bicValue1 = 1e20;  lambdabic1 = c();

  a_vector    = c()
  b_vector    = c()
  sumq_vector = c()

  if(is.null(bic_factor)){
    bic_factor  = 0.3
  }

  all_est = list()

  for( num_pca in num_pca_vec){
    # FPCA
    tuPara[6]   = num_pca
    x_recv_pca  = fda::pca.fd(x_recv, nharm = num_pca)
    x_pca_score = x_recv_pca$scores
    Apena       = create_penamatrix(n_b, num_pca)
    A           = cbind(spMatrix(nrow(Apena$A), ncol(Z)), Apena$A)

    # initial value, estimate subject by subject.
    est_initial = initial_pena(data_list, num_pca, A, est_fix_eff=est_fix_eff, est_treat = est_treat)
    theta       = est_initial$theta
    eta         = est_initial$eta
    lagU        = est_initial$lagU

    #cat(sprintf("num_pca %d length(theta) %d length(eta) %d length(lagU) %d dim(A) (%d %d) \n", num_pca, length(theta), length(eta), length(lagU), dim(A)[1], dim(A)[2]))

    for(lam_g in lamGroup){
      tuPara[1] = lam_g
      bicRes    = fdagrouping(y, Z, x_pca_score, struT, A, theta, eta, lagU, tuPara, precision)

      if(any(is.na(bicRes$beta)) | any(is.nan(bicRes$beta)) | any(is.infinite(bicRes$beta))){
        next
      }
      # spatialgroup::plot_bicres(bicRes)
      theta = bicRes$beta
      eta   = bicRes$eta
      lagU  = bicRes$lagU

      # effective number of parameters
      tempgroup  = adjust_group(group_ind(n_b, num_pca, eta, Apena$Delta, thresh = 1e-4))
      sumq = length(unique(round(theta, 2)))
      sumq = ncol(Z) + length(unique(tempgroup)) * num_pca

      tbicValue   = bicRes$aa + bic_factor *log(n) * (sumq)/n  #*log(log(n+p))

      a_vector    = c(a_vector, bicRes$aa)
      b_vector    = c(b_vector, log(n) * (sumq)/n)
      sumq_vector = c(sumq_vector, sumq)


      if(tbicValue[1] < eval(parse(text = paste("bicValue", 1, sep = "")))){
        assign(paste("est", 1, sep = ""), pack_pena(n_b, ncol(z), x_recv_pca, theta, tempgroup, est_fix_eff=est_fix_eff, est_treat = est_treat))
        assign(paste("bicValue", 1, sep = ""), tbicValue[1])
        assign(paste("lambdabic", 1, sep = ""), c(num_pca, lam_g))
      }

      # store all the results
      all_est[[length(all_est)+1]] = pack_pena(n_b, ncol(z), x_recv_pca, theta, tempgroup, est_fix_eff=est_fix_eff, est_treat = est_treat)

    }
  }
  return(list(est1 = est1, bicValue1 = bicValue1, lambdabic1 = lambdabic1,
              a_vector = a_vector, b_vector = b_vector, all_est = all_est))
}




#' Estimate a partially functional linear regression model with latent group structures using penalization methods.
#' 
#' The estimator is a variant of the penalization estimator implemented with a ADMM algorithm. 
#' The number of groups and the number of principal components are determined by the BIC criterion.
#' 
#' @param data_list A list of data. Several elements must be present in the list. The reponse \code{y},
#' the functional covariate \code{x_recv}, the scalar covariates \code{z}, and an index matrix \code{index}. 
#' The functional covariate \code{x_recv} must be generated from the \code{fda} package by, e.g., spline smoothing. 
#' The scalar covariates \code{z} is a matrix. The index matrix \code{index} is a data.frame recording the structure of the
#' data. The first column of \code{index} is the family number, the second column is the within family index. The column names of \code{index}
#' must be \code{ind_b} and \code{ind_w}.
#' 
#' @param num_pca A vector of candidate number of principal components.
#' 
#' @param penalty A string indicates which penalty function to use, SCAD or LASSO.
#' 
#' @param est_fix_eff A logical value. If \code{TRUE}, then the fixed effects are estimated. Otherwise,
#' the fixed effects are not estimated
#' 
#' @param precision A vector of thresholds for terminating the ADMM algorithm.
#' 
#' @param lam_max_interval A vector indicating the search limits for the maximum of lambda such all the families are clutered into one group.
#' 
#' @param lam_len The size of grid for searching lambda.
#' 
#' @param lamGroup_a A vecotr of candidate values for lambda. If not provided, the algorithm will determine it automatically.
#' 
#' @param scale_pena A logical value, whether to scale the covariates before applying the ADMM algorithm.
#' 
#' @export
pena_est_fda_scale <- function(data_list, num_pca_vec, penalty = 'SCAD', est_fix_eff=TRUE, precision = c(1e-5, 1e-5),
                              lam_max_interval = c(0.1, 10), lam_len = 100, lamGroup_a = NULL, scale_pena = FALSE){
  # this function scales the principal scores
  y      = data_list$y
  #x_dis  = data_list$x_dis
  x_recv = data_list$x_recv
  #t_grid = data_list$t_grid
  z      = data_list$z
  index  = data_list$index
  n_b    = length(unique(index[,1]))
  struT  = as.numeric((index %>% dplyr::group_by(ind_b) %>% dplyr::summarise(n_w=n()))$n_w)
  n      = nrow(data_list$index) # total number of observations

  est_treat = FALSE
  
  tuPara[1]=2        # lambda for group
  tuPara[2]=1        # v
  tuPara[3]=3.0      # gamma parameter in the penalty function
  if(penalty == 'SCAD'){
    tuPara[4]=2        # 1 for LASSO, 2 for SCAD,  3 or others for MCP, for variable selction and grouping
  } else {
    tuPara[4]=2
  }
  tuPara[5]=100000   # maximum iteration numbers
  tuPara[6]=2        # number of pca
  
  
  #cat("In the function pena_est_fda\n")
  #saveRDS(data_list, file='fail.rds')

  # create the pseudo design matrix
  des_dummy   = create_design(index, est_fix_eff=est_fix_eff, mu=TRUE, est_treat = FALSE) # pseudo design, include mean, fixed effects and treatment effects
  Z           = cbind(des_dummy, z)

  # store the estimators
  est1 = c();    bicValue1 = 1e20;  lambdabic1 = c();

  a_vector    = c()
  b_vector    = c()
  sumq_vector = c()
  g_matrix    = c()

  
    bic_factor  = 1
  


  for( num_pca in num_pca_vec){
    # FPCA
    tuPara[6]   = num_pca
    x_recv_pca  = fda::pca.fd(x_recv, nharm = num_pca)
    x_pca_score = x_recv_pca$scores

    # scale the pca score subject-wise?
    if(scale_pena){
    scale_temp = scale(x_pca_score, center = FALSE)
    scale_vec  = attr(scale_temp, 'scaled:scale')
    x_pca_score_scale = scale_temp
    } else {
      x_pca_score_scale = x_pca_score
      scale_vec         = rep(1, ncol(x_pca_score))
    }

    Apena       = create_penamatrix(n_b, num_pca)
    A           = cbind(spMatrix(nrow(Apena$A), ncol(Z)), Apena$Delta %x% Diagonal(x = 1/scale_vec)) # add some zero to account for the finite components

    # initial value, estimate subject by subject.
    est_initial = initial_pena(data_list, num_pca, A, est_fix_eff=est_fix_eff, est_treat = FALSE)
    theta       = est_initial$theta
    theta[(ncol(Z)+1):length(theta)] = theta[(ncol(Z)+1):length(theta)] / rep(scale_vec, n_b)
    eta         = est_initial$eta
    lagU        = est_initial$lagU

    #cat(sprintf("num_pca %d length(theta) %d length(eta) %d length(lagU) %d dim(A) (%d %d) \n", num_pca, length(theta), length(eta), length(lagU), dim(A)[1], dim(A)[2]))

    if(is.null(lamGroup_a)){ # if lamGroup is not available, then find the maximum value of lambda by searching.
      lam_max = find_max_lam(lam_max_interval, n_b, num_pca_vec[length(num_pca_vec)], Apena$Delta, y, Z, x_pca_score_scale,
                             struT, A, theta, eta, lagU, tuPara, precision)
      lamGroup = seq(lam_max + 0.5, 0.1, len = lam_len)
    } else {
      lamGroup = lamGroup_a
    }

    for(lam_g in lamGroup){
      tuPara[1] = lam_g
      bicRes    = fdagrouping(y, Z, x_pca_score_scale, struT, A, theta, eta, lagU, tuPara, precision)

      if(any(is.na(bicRes$beta)) | any(is.nan(bicRes$beta)) | any(is.infinite(bicRes$beta))){
        next
      }
      # spatialgroup::plot_bicres(bicRes)
      theta = bicRes$beta
      theta[(ncol(Z)+1):length(theta)] = theta[(ncol(Z)+1):length(theta)] / rep(scale_vec, n_b)
      eta   = bicRes$eta
      lagU  = bicRes$lagU

      # effective number of parameters
      tempgroup  = adjust_group(group_ind(n_b, num_pca, eta, Apena$Delta, thresh = 1e-1))
      sumq = length(unique(round(theta, 2)))
      sumq = ncol(z) + length(unique(tempgroup)) * num_pca  + as.numeric(est_treat) + as.numeric(est_fix_eff) * n_b

      tbicValue   = bicRes$aa + bic_factor *log(n) * (sumq)/n  # various BIC
      tbicValue[1] =  bicRes$aa + 2 * (sumq)/n               # AIC

      a_vector    = c(a_vector, bicRes$aa)
      b_vector    = c(b_vector, log(n) * (sumq)/n)
      sumq_vector = c(sumq_vector, sumq)
      g_matrix    = rbind(g_matrix, tempgroup)

      for(i in 1:length(bic_factor)){
        if(tbicValue[i] < eval(parse(text = paste("bicValue", i, sep = "")))){
          assign(paste("est", i, sep = ""), pack_pena(n_b, ncol(z), x_recv_pca, theta, tempgroup, est_fix_eff=est_fix_eff, est_treat = est_treat))
          assign(paste("bicValue", i, sep = ""), tbicValue[i])
          assign(paste("lambdabic", i, sep = ""), c(num_pca, lam_g))
        }
      }
    }
  }
  return(list(est1 = est1, bicValue1 = bicValue1, lambdabic1 = lambdabic1))
}

#' Estimate a partially functional linear regression model with latent group structures using a two-step method.
#' 
#' In the first step, a family-wise estimator is calculated. In the second step, the classic K-means algorithm
#' is applied to cluster the families.
#' The number of groups and the number of principal components are determined by the BIC criterion.
#' 
#' @param data_list A list of data. Several elements must be present in the list. The reponse \code{y},
#' the functional covariate \code{x_recv}, the scalar covariates \code{z}, and an index matrix \code{index}. 
#' The functional covariate \code{x_recv} must be generated from the \code{fda} package by, e.g., spline smoothing. 
#' The scalar covariates \code{z} is a matrix. The index matrix \code{index} is a data.frame recording the structure of the
#' data. The first column of \code{index} is the family number, the second column is the within family index. The column names of \code{index}
#' must be \code{ind_b} and \code{ind_w}.
#' 
#' @param num_group_vec A vector of candidate number of groups.
#' 
#' @param num_pca_vec A vector of candidate number of principal components.
#' 
#' @param est_fix_eff A logical value. If \code{TRUE}, then the fixed effects are estimated. Otherwise,
#' the fixed effects are not estimated
#' 
#' 
#' 
#' @export
naive_kmeans <- function(data_list, num_group_vec, num_pca_vec, est_fix_eff=TRUE){
  # estimate the functional coefficients subject-wise, then apply the kmeans algorithm
  # select the number of groups by the BIC criterion.

  x_recv = data_list$x_recv
  n_b    = length(unique(data_list$index[,1]))
  n      = nrow(data_list$index)
  
    bic_factor  = 1
    est_treat = FALSE
  

  # estimate subject-wisely, note that the number of PCs is selected in this step.
  est_subje = family_wise(data_list, num_pca_vec, est_fix_eff=est_fix_eff)
  num_pca   = est_subje$lambdabic1

  # kmeans
  sub_est = est_subje$est1
  xhat_subjectwise = sub_est$xhat_group

  # we need the principal functions to recover the functional coefficients.
  x_recv_pca  = fda::pca.fd(x_recv, nharm = num_pca)

  est1 = NULL; bicValue1 = 1e20; lambdabic1 = 1

  a_vector    = c()
  sumq_vector = c()
  b_vector    = c()

  for(num_group in num_group_vec){
    sub_group_t   = kmeans(t(xhat_subjectwise), centers = num_group, nstart = 100)
    group_index_t = adjust_group(sub_group_t$cluster)
    xhat_group_t  = t(sub_group_t$centers)

    # re-estimate with the estimated group structure
    # we estimate using the following function iterate once (max_iter=1) with a given group structure,
    re_est = fdagroup_est_scale(data_list, num_group = num_group, num_pca = num_pca, est_fix_eff=est_fix_eff,
                                   est_treat=FALSE, max_iter = 1, group_index = group_index_t)

    # BIC
    # xhat_subject_t = matrix(0.0, nrow = num_pca, ncol = n_b)
    # uni_group_index = unique(group_index_t)
    # for(ix in 1:n_b){
    #   xhat_subject_t[,ix] = xhat_group_t[, which(group_index_t[ix] == uni_group_index)]
    # }
    #
    # xhat_subject_fd_t = fd(coef=x_recv_pca$harmonics$coefs %*% xhat_subject_t, basisobj = x_recv_pca$harmonics$basis)
    #
    # est_pack_t = list(muhat = sub_est$muhat, fixeffhat = sub_est$fixeffhat, treathat = sub_est$treathat, zhat = sub_est$zhat,
    #      group=group_index_t, xhat_subject = xhat_subject_t, xhat_sub_fd=xhat_subject_fd_t)

    # bic value
    aa   = mean(predict_fdagroup(re_est, data_list)$pred_error^2, na.rm = TRUE) # RSS
    sumq = num_pca * num_group + ncol(data_list$z) + as.numeric(est_fix_eff) * n_b + as.numeric(est_treat) # number of parameters.

    #bic_factor  = seq(from = 0.30, 1.0, length.out = 8)
    tbicValue   = aa + bic_factor *log(n) * (sumq)/n  #*log(log(n+p))
    tbicValue[1] =  aa + 2 * (sumq)/n               # AIC

    a_vector    = c(a_vector, aa)
    b_vector    = c(b_vector, log(n) * (sumq)/n)
    sumq_vector = c(sumq_vector, sumq)

    for(i in 1:length(bic_factor)){
      if(tbicValue[i] < eval(parse(text = paste("bicValue", i, sep = "")))){
        assign(paste("est", i, sep = ""), re_est)
        assign(paste("bicValue", i, sep = ""), tbicValue[i])
        assign(paste("lambdabic", i, sep = ""), c(num_pca, num_group))
      }
    }
  }
  return(list(est1 = est1, bicValue1 = bicValue1, lambdabic1 = lambdabic1))
}


