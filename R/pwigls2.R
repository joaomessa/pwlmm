#  Fit linear mixed-effects models incorporating sampling weights
#
#  Uses: probability-weighted IGLS (PWIGLS)
#' @useDynLib pwlmm
#'
.onUnload <- function (libpath) {library.dynam.unload("pwlmm", libpath)} #Why it works (using to avoid error message)

#' @export
print.univariate <- function(x, ...){
  cat(noquote("Call:\n"))
  print(x$call)
  cat("\n")
  cat(noquote("Number of observations:"), x$num_obs_1, "\n")
  cat(noquote("Number of groups:"), x$num_obs_2, "\n")
  cat("\n")
  cat(noquote("Fixed effects:\n"))
  result_beta=data.frame(estimate = x$beta$coefficients, std_error = x$beta$standard_errors, t_value = x$beta$coefficients/x$beta$standard_errors)
  names(result_beta) <- c("Estimate", "Std. Error", "t value")
  print(result_beta, digits=4)

  cat("\n")
  cat(noquote("Random effects:\n"))
  cat("\n")
  cat(noquote("Variance components:\n"))
  result_theta = data.frame(estimate = c(x$theta$level2_variances,x$theta$level1_variance) ,
                            se = c(x$theta$se_level2_variances, x$theta$se_level1_variance))
  names(result_theta) <- c("Estimate","Std. Error")
  print(result_theta, digits = 4)

  if(!is.null(x$theta$level2_covariances)){
    var_cov <- x$theta$var_cov
    ujs=nrow(var_cov)
    var_cov <- as.data.frame(format(var_cov, digits=4))
    var_cov <- as.matrix(var_cov)
    var_cov[upper.tri(var_cov, diag = T)] <- ""
    var_cov <- as.data.frame(var_cov)
    var_cov <- var_cov[-1,]
    names(var_cov)[length(names(var_cov))]<-""
    cat("\n")
    cat(if (ujs == 2) noquote("Residuals covariance:\n") else noquote("Residuals covariances:\n"))
    print(var_cov, digits = 4)}

  cat("\n")
  cat(noquote("Note: robust standard errors"))
}

#' @export
residuals.univariate <- function(object, ...){
  object$individual_residuals
}

#' @export
fitted.univariate <- function(object, ...){
  object$fitted_values
}

#' Fit Weighted Linear Multilevel Model
#'
#' Fit a probability-weighted two-level linear model with unequal selection probabilities at each level, via IGLS algorithm.
#'
#' Follows estimation process described in Pfeffermann et al. (1998). Uses probability-weighted IGLS with scaled weights.
#'
#' @param formula a two-sided linear formula object describing both the fixed-effects and random-effects part of the model, with the response on the left of a ~ operator and the terms, separated by + operators, on the right. Random-effects terms are distinguished by vertical bars (|) separating expressions for design matrices from grouping factors.
#' @param data an optional data frame containing the variables in \code{formula}. If not found in data, the variables are taken from the environment of \code{formula} (if specified as a formula) or from the parent frame (if specified as a character vector).
#' @param wj a vector of sampling weights for level two units. Level two units are selected with inclusion probabilities. Then, sampling weights for the level two units are defined as the inverse of these probabilities.
#' @param wij a vector of sampling weights for level one units. After selecting a level two unit, level one units belonging to them are selected with inclusion probabilities. Then, sampling weights for the level one units are defined as the inverse of these probabilities.
#'
#' @references
#' D. Pfeffermann; C. J. Skinner; D. J. Holmes; H. Goldstein; J. Rasbash, 2008,
#' Weighting for Unequal Selection Probabilities in Multilevel Models
#' Journal of the Royal Statistical Society. Series B (Statistical Methodology),
#' Vol. 60, No. 1. (1998), pp. 23-40.
#'
#' @examples
#' pwigls2( Y ~ X1 + X2 + (1|PSU), data = exampledata, wj, wij)
#'
#' @importFrom stats model.matrix model.response
#' @return Estimated list of estimators
#' @export
pwigls2 <-function(formula, data = NULL, wj, wij){
  clprint <- cl <- match.call()
  ma <- match(c("formula", "data", "wj", "wij"), names(cl), 0L)
  cl <- cl[c(1L, ma)]
  cl$drop.unused.levels <- TRUE
  cl[[1L]] <- quote(stats::model.frame)
  cl$formula <- lme4::subbars(formula)
  fr <- eval(cl, parent.frame())
  fr <- lme4::factorize(cl$formula, fr, char.only = TRUE)
  attr(fr, "formula") <- formula

  vars_model <- formula[[length(formula)]]
  z <- lme4::findbars(vars_model)[[1]]
  fr <- fr[order(fr[,deparse(z[[3]])]),]

  x <- model.matrix(eval(substitute( ~ foo, list(foo = vars_model[[2]]))), fr)
  clusters<-as.matrix(cbind(fr[,deparse(z[[3]])]))
  z <- model.matrix(eval(substitute( ~ foo, list(foo = z[[2]]))), fr)
  m <- nrow(as.matrix(unique(clusters)))

  y <- model.response(fr)
  n <- length(y)

  wj<-fr$`(wj)`
  wi_j<-fr$`(wij)`

  name1 <- colnames(x)
  name2 <- colnames(z)
  q <- ncol(z)
  s <- ((q*(q+1))/2)+1
  one_to_sminus1 <- 1:(s-1)
  matrix_invvech <- invvech_eigen(one_to_sminus1)

  if(q != 1){
    H <- t(diag(s)[matrix_invvech,])
  } else {
    H <- cbind(diag(s)[matrix_invvech,])
  }

  p <- ncol(x)

  name3 <- matrix(paste ("sigmau_", rep((0:(q-1)), (q:1)), ((sequence(q:1))-1)+ (rep((0:(q-1)), (q:1))), sep=""))
  name3<-rbind(name3, "sigma2_e")

  nepg <- tapply(clusters, clusters, function(x) NROW(x))
  panelsetup <- as.matrix(cbind(nepg, cumsum(nepg)))
  panelsetup[,2] <- panelsetup[,2] - panelsetup[,1]

  #-------------- Calculating the Scaled Weights --------------
  Scaled_W <- scaled_weight(n,
                             m,
                             panelsetup,
                             wi_j,
                             wj)

  TJS <- tjs_uni_beta(p,
                      m,
                      panelsetup,
                      x,
                      y,
                      z,
                      Scaled_W$wi_j_star,
                      q,
                      s,
                      H,
                      Scaled_W$wj_star)

  beta0 <- solve(TJS$somat1) %*% TJS$somat3
  beta0 <- as.numeric(beta0)

  i_t <- initial_theta(
    beta0,
    x,
    y,
    m,
    panelsetup,
    Scaled_W$wi_j_star,
    Scaled_W$wj_star,
    p)

  sit <- diag(2, q)
  teta0 <- c(sit[lower.tri(sit, diag = T)], i_t$wj_t6/i_t$wj_aux)

  #--------------------------------------------------------------------------
  #                     IGLS - ITERATIVE
  #--------------------------------------------------------------------------

  itera = 0

  ## Trick
  beta_ant = beta0
  beta = beta_ant*2
  teta_ant=teta0
  teta = teta_ant*2

  while (itera<= 200  & (any(abs((teta-teta_ant))> 0.000001) | any(abs((beta-beta_ant))> 0.000001))){

    if (itera == 0) {
      teta <- teta0
    }

    tsit <- teta[s]*sit

    objeto1 <- iterative_uni_beta(p,
                                  m,
                                  tsit,
                                  TJS$T1,
                                  TJS$T2,
                                  TJS$T3,
                                  TJS$T4,
                                  TJS$T5,
                                  Scaled_W$wj_star)

    #----------- beta--------------

    if (itera != 0){
      beta_ant = beta
    }

    solve_s_matp <- solve(objeto1$s_matp)
    beta = solve_s_matp %*% objeto1$s_matq
    beta <- as.numeric(beta)

    objeto2 <- iterative_uni_theta(s,
                                   p,
                                   m,
                                   sit,
                                   tsit,
                                   panelsetup,
                                   y,
                                   x,
                                   z,
                                   q,
                                   beta,
                                   Scaled_W$wj_star,
                                   teta[s],
                                   objeto1$AJS,
                                   TJS$T5,
                                   Scaled_W$wi_j_star,
                                   TJS$H_K,
                                   TJS$TR_T5_HK,
                                   TJS$T5_HK)

    #-----theta ----

    if (itera != 0) {
      teta_ant = teta
    }

    solve_r_mat <- solve(objeto2$r_mat)
    teta = solve_r_mat %*% objeto2$s_mat
    teta <- as.numeric(teta)

    ie_teta <- invvech_eigen(teta[one_to_sminus1]) #Variance covariance matrix sigma_u
    sit <- solve(ie_teta)

    #------End of iterative process----------
    itera = itera + 1
  }
  v=diag(ie_teta)
  teta_var1 <- teta[s]
  if(any(v<0) | teta_var1<0)
    stop('Model failed to converge: negative variance component(s)')

  #-------------------------------------------------------------------------

  #-------------------------------------------------------------------------
  # Variances
  #-------------------------------------------------------------------------

  variances_residuals <- uni_variances_residuals(
    p,
    s,
    beta,
    teta,
    y,
    x,
    m,
    panelsetup,
    TJS$T2,
    TJS$T5,
    Scaled_W$wj_star,
    sit,
    n,
    q,
    z,
    ie_teta,
    teta_var1,
    Scaled_W$wi_j_star,
    TJS$H_K,
    TJS$TR_T5_HK,
    TJS$T5_HK)

  var_beta = solve_s_matp%*%((m /(m-1))*variances_residuals$s_matc)%*%solve_s_matp
  dp_beta = sqrt(diag(var_beta))

  #var_teta = 2*solve(r_mat) in other context
  var_teta = solve_r_mat%*%(m/(m-1)*variances_residuals$s_matd)%*% solve_r_mat
  dp_teta = sqrt(diag(var_teta))

  #/*-------------------------------------------------------------------------
  #  Residuals
  #-------------------------------------------------------------------------*/

  if(any(variances_residuals$var_u<0))
    warning('Model failed to converge: negative residual(s) variance(s)')

  dp_u <- sqrt(variances_residuals$var_u)

  colnames(variances_residuals$u) <- colnames(dp_u) <-paste("u", 0:(q-1), sep="")
  rownames(variances_residuals$u) <- rownames(variances_residuals$var_u) <- as.vector(unique(clusters))
  names(variances_residuals$v) <- names(variances_residuals$yhat) <- rownames(fr)
  names(teta) <- name3

  diag_matrix_invvech <- diag(matrix_invvech)
  teta_var2 <- teta[diag_matrix_invvech]
  se_teta_var2 <- dp_teta[diag_matrix_invvech]

  se_teta_var1 <- dp_teta[s]

  names(teta_var1) <- "Residual"
  names(teta_var2) <- colnames(ie_teta) <- rownames(ie_teta) <- name2

  if(s!=2){
    minus_dmi <- (one_to_sminus1)[-diag_matrix_invvech]
    teta_cov <- teta[minus_dmi]
    se_cov <- dp_teta[minus_dmi]
  } else {
    teta_cov <- se_cov <- NULL
  }

  names(beta) <- name1
  list_igls <-list(beta = list(coefficients = beta, standard_errors = dp_beta),
                   theta = list(level2_variances = teta_var2, se_level2_variances = se_teta_var2, level2_covariances = teta_cov,
                                se_level2_covariances = se_cov, var_cov = ie_teta, level1_variance = teta_var1,
                                se_level1_variance = se_teta_var1),
                   num_obs_1 = n, num_obs_2 = m, call = clprint, fitted_values = variances_residuals$yhat,
                   individual_residuals = variances_residuals$v,
                   group_residuals = list(coefficients = variances_residuals$u, standard_errors = dp_u),
                   iterations = itera
  )
  class(list_igls) <- 'univariate'
  list_igls
}

#' @export
print.multivariate <- function(x, ...){
  cat(noquote("Call:\n"))
  print(x$call)
  cat("\n")
  cat(noquote("Number of observations:"), x$num_obs_1, "\n")
  cat(noquote("Number of groups:"), x$num_obs_2, "\n")
  cat(noquote("Number of repeated measures:"), x$num_time_obs, "\n")
  cat("\n")
  cat(noquote("Fixed effects:\n"))
  result_beta=data.frame(estimate = x$beta$coefficients, std_error = x$beta$standard_errors, t_value = x$beta$coefficients/x$beta$standard_errors)
  names(result_beta) <- c("Estimate", "Std. Error", "t value")
  print(result_beta, digits = 4)

  cat("\n")
  cat(noquote("Random effects:\n"))
  result_theta=data.frame(estimate = x$theta$coefficients, std_error = x$theta$standard_errors)
  names(result_theta) <- c("Estimate", "Std. Error")
  print(result_theta, digits = 4)

  cat("\n")
  cat(noquote("Note: robust standard errors"))
}

#' @export
residuals.multivariate <- function(object, ...){
  object$individual_residuals
}

#' @export
fitted.multivariate <- function(object, ...){
  object$fitted_values
}

#' @importFrom stats toeplitz
teta_structure_toep <- function(teta, s){
    return(toeplitz(as.numeric(teta[2:s])))
}
teta_structure_uns <- function(teta, s, tt, itera){
  teta_uns <- matrix(0,tt,tt)
  teta_uns[lower.tri(teta_uns, diag=TRUE)] <- as.numeric(teta[2:s])
  teta_uns <- as.matrix(Matrix::forceSymmetric(teta_uns,uplo="L"))
  if (itera == 1)
    diag(teta_uns) <- as.numeric(teta[2])
  return(teta_uns)
}

delta_gen <- function(rot, s){
  tt <- sum(rot)
  k=0
  lag_size = length(rot)-1
  lag_list=replicate(lag_size,NULL)
  for(i in 1:lag_size)
    if(rot[i]==1){
      k=k+1
      l = k
      for(j in (i+1):(lag_size+1)){
        if(rot[j] == 1){
          l = l+1
          lag_dist = j - i
          lag_list[[lag_dist]] = c(lag_list[[lag_dist]],k,l)
        }
      }
    }

  lag_list=Filter(Negate(is.null), lag_list)

  lag_m=diag(tt)
  delta_matrix2=list(lag_m)
  for(i in seq_along(lag_list)){
    for(j in seq(1,length(lag_list[[i]]),by = 2)){
      lag_m[lag_list[[i]][j],lag_list[[i]][j+1]] <- 1
    }
    delta_matrix2[[i+1]] <- lag_m
    lag_m=diag(tt)
  }

  delta_matrix2 <- lapply(delta_matrix2,function(x) Matrix::forceSymmetric(x,uplo="U"))
  DELTA=c(list(matrix(0,tt,tt)),lapply(delta_matrix2, function(x) as.matrix(x)))
  return(list(delta = DELTA,
              teta_loop_genlin = quote(Reduce('+', mapply("*",DELTA[-1], teta[-1], SIMPLIFY = F))),
              name_type = matrix(paste("Genlin ", 1:(length(DELTA)-1), sep=""))))
}

#There's a more effective way which doesn't make use of 'for'
delta_toep <- function(tt){
  k=0
lag_size=tt-1
lag_list=replicate(lag_size,NULL)
for(i in 1:lag_size){
  k=k+1
  l=0
  for(j in i:(lag_size)){
    l=l+1
    lag_list[[i]] = c(lag_list[[i]],l,l+k)
  }
}
lag_m=matrix(0,tt,tt)
delta_matrix2=list(diag(tt))
for(i in seq_along(lag_list)){
  for(j in seq(1,length(lag_list[[i]]),by = 2)){
    lag_m[lag_list[[i]][j],lag_list[[i]][j+1]] <- 1
  }
  delta_matrix2[[i+1]] <- lag_m
  lag_m=matrix(0,tt,tt)
}

delta_matrix2 <- lapply(delta_matrix2,function(x) Matrix::forceSymmetric(x,uplo="U"))
DELTA=c(list(lag_m),lapply(delta_matrix2, function(x) as.matrix(x)))
return(list(delta = DELTA, teta_loop_toep = quote(teta_structure_toep(teta, s)),
            name_type = matrix(paste("TOEP ", 1:tt, sep=""))))
}

delta_uns <- function(tt){
  k=0
  tamanho_delta <- sum(1:tt)
  delta_matrix2 <- rep(list(matrix(0,tt,tt)), tamanho_delta)
  for(i in 1:tt)
    for(j in i:tt){
      k=k+1
      delta_matrix2[[k]][j,i] <- 1
    }
  delta_matrix2 <- lapply(delta_matrix2,function(x) Matrix::forceSymmetric(x,uplo="L"))
  DELTA=c(list(matrix(0,tt,tt)),lapply(delta_matrix2, function(x) as.matrix(x)))
  return(list(delta = DELTA, teta_loop_uns = quote(teta_structure_uns(teta, s, tt, itera)),
              name_type = matrix(paste("UNS ", 1:sum(1:tt), sep=""))))
}
