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
    K = K_calculate(X_1,X_2, cov = cov)
  }

  observed_test_stats = calc_ky_fan_k(cov_cor_X_1,cov_cor_X_2,K)

  # Create empty matrix of permutation statistics we will fill
  permutation_stat_matrix = matrix(data=0,nrow=n_perm,ncol=length(observed_test_stats))
  colnames(permutation_stat_matrix) = names(observed_test_stats)

  # Permute data and record test statistics using permuted data in permutation_stat_matrix
  X_1_X_2 = cbind(X_1,X_2)

  set.seed(seed)

  for(j in 1:n_perm){
    # Randomly permute the data
    permuted_order = sample(1:ncol(X_1_X_2))
    permuted_group_1 = X_1_X_2[,permuted_order[1:ncol(X_1)]]
    permuted_group_2 = X_1_X_2[,permuted_order[(ncol(X_1)+1):(ncol(X_1)+ncol(X_2))]]

    # Calculate test statistics for permuted data, and add to permutation_stat_matrix
    if(cov == TRUE){
      cov_cor_permuted_group_1 = cov(t(permuted_group_1))
      cov_cor_permuted_group_2 = cov(t(permuted_group_2))
    }else{
      cov_cor_permuted_group_1 = cor(t(permuted_group_1))
      cov_cor_permuted_group_2 = cor(t(permuted_group_2))
    }

    permuted_test_stat = calc_ky_fan_k(cov_cor_permuted_group_1,cov_cor_permuted_group_2,K)
    permutation_stat_matrix[j,] = c(permuted_test_stat)
  }

  if(min_P == TRUE){
    # Calculate p values for statistics in each permutation, and then take minimum p across these statistics to get
    # permutation distribution of minimum p value. Note here we do not need to normalize the statistics since the same
    # normalization is applied to every Ky-Fan(k) norm so this will not affect the p-value.
    permutation_p_value_matrix = (n_perm-apply(permutation_stat_matrix,2,rank,ties.method='min')+1)/n_perm

    # Calculate observed p values
    observed_p_values = (colSums(t(t(permutation_stat_matrix)>=observed_test_stats))+1)/(n_perm+1)

    # Add minimum p value column to permutation value matrix and calculate p-value of minimum p-value from observed data
    permutation_min_p_values = apply(permutation_p_value_matrix,1,min)
    observed_min_p_value = min(observed_p_values)
    RACT_p_value = (sum(permutation_min_p_values <= observed_min_p_value) + 1)/(n_perm + 1)

  }else{
    # Normalize observed test statistics, and permutation test statistics by means and variances estimated via permutation
    test_stat_means = apply(permutation_stat_matrix,2,mean)
    test_stat_sd = apply(permutation_stat_matrix,2,sd)

    normalized_observed_test_stats = (observed_test_stats - test_stat_means)/test_stat_sd
    normalized_permutation_stat_matrix = sweep(permutation_stat_matrix,2,test_stat_means,'-')
    normalized_permutation_stat_matrix = sweep(normalized_permutation_stat_matrix,2,test_stat_sd,'/')

    # Calculate observed p values and RACT p value
    observed_p_values = ((colSums(t(t(normalized_permutation_stat_matrix)>=normalized_observed_test_stats))+1)/(n_perm+1))
    max_normalized_permutation_stats = apply(normalized_permutation_stat_matrix,1,max)
    RACT_p_value = (sum(max_normalized_permutation_stats>=max(normalized_observed_test_stats))+1)/(n_perm+1)
  }

  return(list('RACT p value'=RACT_p_value, 'Individual Ky-Fan(k) p values'= observed_p_values))
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
#' @noRd
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
#' @export
#' @examples
#' set.seed(1)
#' n = 50
#' p = 250
#' a_1 = matrix(rnorm(p),ncol=1)
#' X_1 = a_1%*%t(a_1)%*%matrix(rnorm(n*p),nrow=p,ncol=n)
#' a_2 = matrix(rnorm(p),ncol=1)
#' X_2 = a_2%*%t(a_2)%*%matrix(rnorm(n*p),nrow=p,ncol=n)
#' K_calculate(X_1,X_2)
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
  # Reject
  n = 50
  p = 250

  set.seed(1)
  a_1 = matrix(rnorm(p),ncol=1)
  X_1 = a_1%*%t(a_1)%*%matrix(rnorm(n*p),nrow=p,ncol=n)

  a_2 = matrix(rnorm(p),ncol=1)
  X_2 = a_2%*%t(a_2)%*%matrix(rnorm(n*p),nrow=p,ncol=n)

  K = NULL

  # Small reject
  n = 50
  p = 250
  snr = .02

  set.seed(1)
  X_1 = matrix(rnorm(n*p),nrow=p,ncol=n)

  cov_X_2 = diag(p)*(1-snr)+snr
  svd_cov_X_2 = svd(cov_X_2)
  A = svd_cov_X_2$u%*%diag(sqrt(svd_cov_X_2$d))%*%t(svd_cov_X_2$v)
  X_2 = A%*%matrix(rnorm(n*p),nrow=p,ncol=n)

  K = NULL

}

