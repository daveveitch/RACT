#' RACT method
#'
#' @description This function takes in two data matrices, and tests for the equality of covariance using rank-adaptive covariance testing.
#'
#' @param X_1 p x n data matrix (p dimensions, n observations) for group 1
#' @param X_2 p x n data matrix (p dimensions, n observations) for group 2
#' @param K If NULL, Ky-Fan(k) norms from 1 to K are included where K smallest integer such that sum of top K singular values of the pooled covariance
#'   are greater than 80% of the sum of all of the singular values. Otherwise this argument takes a vector of integers, and these integers represent the
#'   Ky-Fan(k) norms that are included. E.g. K=c(1,5) means RACT is based on the Ky-Fan(1) and Ky-Fan(5) norm.
#' @param min_P Whether a minimum-P value statistic is used.
#' @param cov If cov=TRUE RACT tests for the equality of the covariance of X_1,X_2. If cov=FALSE RACT tests for the equality of the correlation of X_1,X_2.
#' @param seed Allows user to set a seed for the permutation procedure.
#' @return a list of length 2. 'adaptive_p_value' is the p-value associated with RACT's adaptive test statistic. 'p_value_vector' a vector with p-values
#'   associated with individual Ky-Fan(k) norms
#' @export
#'
#' @examples
#' X_1 =
RACT <- function(X_1,X_2,n_perm=1000,K=NULL,min_P=FALSE,cov=TRUE,seed=NULL){

  # Calculate covariance/correlation matrices for each group
  if(cov == TRUE){
    cov_cor_X_1 = cov(t(X_1))
    cov_cor_X_2 = cov(t(X_2))
  }else if(cov == FALSE){
    cov_cor_X_1 = cor(t(X_1))
    cov_cor_X_2 = cor(t(X_2))
  }else{
    stop('cov must be set to TRUE/FALSE')
  }

  if(is.null(K)){

  }

  observed_test_stats = calc_ky_fan_k(cov_cor_X_1,cov_cor_X_2,K)

  permutation_stat_matrix = calc_permutation_stat_matrix(cbind(X_1,X_2),
                                                         ncol(X_1),
                                                         ncol(X_2),
                                                         K,n_perm,
                                                         observed_test_stats,
                                                         cov_cor=cov_cor)

  permutation_p_value_matrix = calc_permutation_p_value_matrix(permutation_stat_matrix)
  permutation_p_value_matrix = permutation_add_min_p_cols(permutation_p_value_matrix)

  observed_p_values = calc_observed_p_values(observed_test_stats,permutation_stat_matrix)
  observed_p_values = observed_add_min_p_cols(observed_p_values,
                                              permutation_p_value_matrix)

  reject_matrix = (observed_p_values <= alpha)*1

  return(reject_matrix)
}

#' Calculate Ky-Fan(k) norm
#'
#' @description This function takes in two covariance/correlation matrices, and calculates the Ky-Fan(k) norm of their differences for k=1,...,K
#'
#' @param cov_cor_X_1 p x p (p dimensions) covariance or correlation matrix for group 1
#' @param cov_cor_X_2 p x p  (p dimensions) covariance or correlation matrix for group 2
#' @param K Vector of integers, and these integers represent the Ky-Fan(k) norms that are included.
#'   E.g. K=c(1,5) means RACT is based on the Ky-Fan(1) and Ky-Fan(5) norm.
#' @return named vector of Ky-Fan(k) norms where names of each entry are values of k
calc_ky_fan_k<-function(cov_cor_X_1,cov_cor_X_2,K){
  diff_mat = cov_cor_X_1-cov_cor_X_2

  # Calculate top K singular values of difference matrix
  top_K_singular_values = sort(RSpectra::svds(diff_mat,k=max(K),nu=0,nv=0)$d,decreasing=TRUE)

  ky_fan_k_norms = cumsum(top_K_singular_values)[K]
  names(ky_fan_k_norms) = K

  return(ky_fan_k_norms)
}

#' This function calculates the smallest K, such that the sum of the top K singular values of the covariance/correlation
#'   matrix of the data is at least K_pct of the sum of all singular values.
#'
#' @description This function takes in two covariance/correlation matrices, and calculates the Ky-Fan(k) norm of their differences for k=1,...,K
#'
#' @param X_1 p x n data matrix (p dimensions, n observations) for group 1
#' @param X_2 p x n data matrix (p dimensions, n observations) for group 2
#' @param K_pct 0 < K_pct <= 1, where K_pct determines how large K should be so that sum of top K singular values of pooled covariance/correlation at least
#'   K_pct of the sum of all singular values. Default set at 0.8.
#' @param cov If cov=TRUE then K based off of pooled covariance, if cov=FALSE then based off of pooled correlation.
#'
#' @return vector of values c(1,...,K)
K_calculate<-function(X_1,X_2,K_pct = 0.8,cov = TRUE){

  combined_X = cbind(X_1,X_2)

  if(cov == TRUE){
    cov_cor_combined_X = cov(t(combined_X))
  }else if(cov == FALSE){
    cov_cor_combined_X = cor(t(combined_X))
  }else{
    stop('cov must be set to TRUE/FALSE')
  }

  # At most pooled correlation/covariance will have n (number of observations) non-zero singular values, avoids situation where this
  # function attempts to take a SVD of an extremely large matrix
  singular_values = sort(RSpectra::svds(cov_cor_combined_X,
                                        k=min(nrow(combined_X),ncol(combined_X)),nu=0,nv=0)$d,decreasing=TRUE)

  # Calculate sum of all singular values and determine where the cutoff should be
  singular_value_cutoff = K_pct * sum(singular_values)
  K = seq(1,which(cumsum(singular_values) >= singular_value_cutoff)[1])

  return(K)
}


# For testing
test_fun<-function(){
  n = 50
  p = 250
  cov_matrices = generate_covariance_matrix('(LowRank)',1,p,5)
  X_1 = x_matrix=t(MASS::mvrnorm(n=n,mu=rep(0,p),Sigm=cov_matrices[[1]]))
  X_2 = x_matrix=t(MASS::mvrnorm(n=n,mu=rep(0,p),Sigm=cov_matrices[[2]]))
  K = NULL
}

