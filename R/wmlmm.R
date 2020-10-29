#  Fit linear mixed-effects models incorporating sampling weights
#
#  Uses: probability-weighted IGLS (PWIGLS)

#' Fit Weighted Multivariate Linear Multilevel Model to Longitudinal Data
#'
#' Fit a two-level probability-weighted multivariate linear model with a linear error covariance matrix structure, via IGLS algorithm.
#'
#' Follows estimation process described in Veiga et al. (2014). Uses probability-weighted IGLS with scaled weights.
#'
#' @param formula a linear formula object with the response on the left of a ~ operator and the terms, separated by + operators, on the right.
#' @param data an optional data frame containing the variables in \code{formula}. If not found in data, the variables are taken from the environment of \code{formula} (if specified as a formula) or from the parent frame (if specified as a character vector).
#' @param ID3 vector of indexes for level two units
#' @param ID2 vector of indexes for level one units.
#' @param ID1 vector of successive measurements within the same level one unit, for all units.
#' @param wj a vector of sampling weights for level two units. Level two units are selected with inclusion probabilities. Then, sampling weights for the level two units are defined as the inverse of these probabilities.
#' @param wij a vector of sampling weights for level one units. After selecting a level two unit, level one units belonging to them are selected with inclusion probabilities. Then, sampling weights for the level one units are defined as the inverse of these probabilities.
#' @param type type of structure imposed in the error covariance matrix; "toep" refers to the toeplitz, "uns" refers to the unestructured and "genlin" refers to the general linear.
#' @param rot vector of 0's and 1's related to the measurements in time when \code{"genlin"} is passed to the \code{type} argument. Use 1 if the data were collected in that specific time unit, and 0 otherwise.
#'
#' @references
#' A. Veiga, P. W. F. Smith and J. J. Brown (2014),
#' The use of sample weights in multivariate multilevel models with an application to income data collected by using a rotating panel survey
#' Journal of the Royal Statistical Society. Series C (Applied Statistics)
#' Vol. 63, No. 1 (JANUARY 2014), pp. 65-84 (20 pages)
#'
#' @examples
#'  wmlmm ( Y ~ X1 + X2, data = exampledata, idd, wave, wj, wi_j, "toep")
#'
#' @importFrom stats model.matrix model.response
#' @return Estimated list of estimators
#' @export
wmlmm <- function(formula, data = NULL, ID3, ID2, ID1, wj, wij, type, rot = NULL){
  clprint <- cl <- match.call()
  ma <- match(c("formula", "data", "ID3", "ID2", "ID1",
                "wj", "wij"), names(cl), 0L)
  name3_1<-match("ID3", names(cl), 0L)
  name3_1<- cl[name3_1]

  cl <- cl[c(1L, ma)]
  cl$drop.unused.levels <- TRUE
  cl[[1L]] <- quote(stats::model.frame)
  cl <- eval(cl, parent.frame())

  cl <- cl[order(cl$`(ID3)`, cl$`(ID2)`, cl$`(ID1)`),]

  mm=model.matrix(formula, cl)

  x <- mm[,-1, drop=FALSE]
  x_names <- colnames(x)
  y <- as.numeric(model.response(cl))

  n <- length(y)
  z <- rep(1, n) #1

  clusters<-as.matrix(cbind(cl$`(ID3)`))
  m <- nrow(as.matrix(unique(clusters)))

  wj <- cl$`(wj)`
  wi_j <- cl$`(wij)`

  if(type=="genlin"){
    tt <- sum(rot)
  } else {
    tt=n/nrow(unique(cbind(cl$`(ID3)`,cl$`(ID2)`)))
  }

  n_nivel_1 = n/tt

  x_waves=rep(list(diag(tt)), n_nivel_1)
  x_waves=do.call(rbind, x_waves)

  x<-as.matrix(cbind(x_waves,x))
  p <- ncol(x)

  if(type=="toep"){
    lista <- delta_toep(tt)
    DELTA <- lista$delta
    teta_loop <- lista$teta_loop_toep
    name_type <- lista$name_type
  } else if(type=="uns") {
    lista <- delta_uns(tt)
    DELTA <- lista$delta
    teta_loop <- lista$teta_loop_uns
    name_type <- lista$name_type
  } else {
    lista <- delta_gen(rot)
    DELTA <- lista$delta
    teta_loop <- lista$teta_loop_genlin
    name_type <- lista$name_type
  }

  s <- length(DELTA)
  h_matrix = diag(s)[1,]

  name1 <- matrix(paste (deparse(substitute(ID1)), 1:tt, sep=" "))
  name1 <- if (is.null(x_names)) name1 else rbind(name1, as.matrix(x_names))
  name3 = rbind(deparse(substitute(ID3)), name_type)

  nepg <- tapply(clusters, clusters, function(x) NROW(x))
  panelsetup <- as.matrix(cbind(nepg, cumsum(nepg)))
  panelsetup[,2] <- panelsetup[,2] - panelsetup[,1]

  #-------------- Calculating the Scaled Weights --------------
  Scaled_W <- scaled_weight(n,
                            m,
                            panelsetup,
                            wi_j,
                            wj)

  inv_wi_j_star <- 1/Scaled_W$wi_j_star

  beta0 <- initial_beta(
    p,
    x,
    y,
    m,
    panelsetup,
    Scaled_W$wi_j_star,
    Scaled_W$wj_star)

  beta0 <- solve(beta0$somat1) %*% beta0$somat3
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

  teta0 = rep(0.5, s)
  teta0[2]=i_t$wj_t6/i_t$wj_aux

  #--------------------------------------------------------------------------
  #                     IGLS - ITERATIVE
  #--------------------------------------------------------------------------

  itera = 1

  ss=s*s

  ## Trick
  beta_ant = beta0
  beta = beta_ant*2
  teta_ant=teta0
  teta=teta_ant * 2

  while (itera<= 200  & (any(abs((teta-teta_ant))> 0.000001) | any(abs((beta-beta_ant))> 0.000001) )){

    if (itera != 1) {
      teta <- teta
    } else{
      teta <- teta0
    }

    teta_matrix <- eval(teta_loop)

    objeto = iterative_multi_beta(
      p,
      m,
      panelsetup,
      x,
      y,
      z,
      Scaled_W$wi_j_star,
      Scaled_W$wj_star,
      as.matrix(teta),
      tt,
      teta_matrix)

    #----------- beta--------------
    solve_s_matp <- solve(objeto$s_matp)

    if (itera != 1){
      beta_ant = beta
    }

    beta <- solve_s_matp %*% objeto$s_matq
    beta <- as.numeric(beta)

    objeto2=iterative_multi_teta(
      m,
      beta,
      panelsetup,
      x,
      y,
      Scaled_W$wi_j_star,
      inv_wi_j_star,
      Scaled_W$wj_star,
      tt,
      DELTA,
      objeto$invvjs,
      h_matrix,
      s,
      ss,
      p)

    matr=colSums(objeto2$R)
    mats=colSums(objeto2$S)

    r_mat = matrix(matr,s)
    s_mat = matrix(mats,s)

    #-----theta = u0j , eij----#

    if (itera != 1) {
      teta_ant = teta
    }

    solve_r_mat <- solve(r_mat)
    teta <- solve_r_mat %*% s_mat
    teta <- as.numeric(teta)

    #------End of iterative process----------
    itera = itera + 1
  }

  teta_matrix <- eval(teta_loop)

  #-------------------------------------------------------------------------
  # Variances
  #-------------------------------------------------------------------------
  teta1 <- teta[1]
  teta_1 <- solve(teta1)

  variances_resids <- multi_variances_residuals(
    p,
    s,
    beta,
    teta,
    y,
    x,
    z,
    Scaled_W$wi_j_star,
    m,
    panelsetup,
    objeto2$R,
    objeto2$S,
    Scaled_W$wj_star,
    tt,
    teta_1,
    teta_matrix,
    n,
    teta1)

  dp_u <- sqrt(variances_resids$var_u)

  var_beta = solve_s_matp%*%((m/(m-1))*variances_resids$s_mat_c)%*%solve_s_matp
  dp_beta = sqrt(diag(var_beta))

  #var_teta = 2*solve(r_mat) in other context
  var_teta = solve_r_mat%*%(m/(m-1)*variances_resids$s_mat_d)%*%solve_r_mat
  dp_teta = sqrt(diag(var_teta))

  Sigma_r=teta1*matrix(1,tt,tt)+ teta_matrix

  names(variances_resids$u) <- names(dp_u) <- as.vector(unique(clusters))
  #names(variances_resids$v) <- names(variances_resids$yhat) <- rownames(mm) necessary? 6x heavier
  names(beta) <- name1
  names(teta) <- name3
  list_genlin <-list(beta = list(coefficients = beta, standard_errors = dp_beta, initial_values = beta0),
                     theta = list(coefficients = teta, standard_errors = dp_teta, initial_values = teta0),
                     mix_teta_matrix = teta_matrix, tot_var_matrix = Sigma_r, num_obs_1 = n_nivel_1, num_obs_2 = m,
                     num_time_obs = tt, call = clprint, fitted_values = variances_resids$yhat,
                     individual_residuals = variances_resids$v,
                     group_residuals = list(coefficients = variances_resids$u, standard_errors = dp_u),
                     iterations = itera-1
  )
  class(list_genlin) <- 'multivariate'
  list_genlin
}
