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
#' @param alpha Prescribed type I error rate of RACT.
#' @return a list of length 2. 'adaptive_p_value' is the p-value associated with RACT's adaptive test statistic. 'p_value_vector' a vector with p-values
#'   associated with individual Ky-Fan(k) norms
#' @export
#'
#' @examples
#' X_1 =
RACT <- function(X_1,X_2,n_perm=1000,K=NULL,min_P=FALSE,cov=TRUE,seed=NULL,alpha=0.05){

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



  observed_test_stats = test_stat_20240925(cov_cor_X_1,cov_cor_X_2,
                                           K,
                                           test_stat_prefixes)

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

calc_ky_fan_k<-function(mat_1,mat_2,max_k,test_stat_prefixes){
  # This function returns the sum of the top 1,2,...,max_k singular values, plus the sum of squares,
  # plus the frobenius norm of the matrix of the difference between the matrices mat_1 and mat_2
  # INPUT
  # mat_1                 - pxp matrix
  # mat_2                 - pxp matrix
  # max_k                 - integer, maximum singular value to compute
  # test_stat_prefixes    - strings indicating what statistics to compute (e.g. sum, sumsq, frob)
  # OUTPUT
  # test_stat_vec - a named vector with calculated test statistics

  diff_mat = mat_1-mat_2

  top_k_singular_values = sort(RSpectra::svds(diff_mat,k=max_k,nu=0,nv=0)$d,decreasing=TRUE)

  test_stat_values = c()
  test_stat_names = c()

  if('sum' %in% test_stat_prefixes){
    test_stat_values = c(test_stat_values, cumsum(top_k_singular_values))
    test_stat_names = c(test_stat_names, paste("sum", 1:max_k, sep="_"))
  }
  if('sumsq' %in% test_stat_prefixes){
    test_stat_values = c(test_stat_values,cumsum(top_k_singular_values**2))
    test_stat_names = c(test_stat_names, paste("sumsq", 1:max_k, sep="_"))
  }
  if('frob' %in% test_stat_prefixes){
    test_stat_values = c(test_stat_values,sum(diff_mat**2))
    test_stat_names = c(test_stat_names, 'frob_all')
  }

  names(test_stat_values) = test_stat_names

  return(test_stat_values)
}


