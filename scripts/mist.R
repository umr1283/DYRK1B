#-------------------------------------------------------------------------------
# Name - MiST: Mixed effects Score Test for continuous outcomes
# Desc - Test for association between a set of SNPS/genes and continuous
#        or binary outcomes by including variant characteristic information
#        and using (weighted) score statistics.
# Original From - https://cran.r-project.org/package=MiST
# References - https://doi.org/10.1002/gepi.21717
#-------------------------------------------------------------------------------

#' mist
#'
#' Test for association between a set of SNPS/genes and continuous outcomes by
#' including variant characteristic information and using score statistics.
#'
#' @param y [numeric] A numeric vector of the continuous outcome variables. Missing values are not allowed.
#' @param X [numeric] A numeric matrix of covariates with rows for individuals and columns for covariates.
#'   If there is no covariate, it does not need to be specified
#' @param G [numeric] A numeric genotype matrix with rows for individuals and columns for SNPs.
#'   Each SNP should be coded as 0, 1, and 2 for AA, Aa, aa, where A is a major allele and a is a minor allele.
#'   Missing genotypes are not allowed.
#' @param Z [numeric]a numeric matrix of second level covariates for variant characteristics.
#'   Each row corresponds to a variant and each column corresponds to a variant characteristic.
#'   If there is no second level covariates, a vector of 1 should be used.
#' @param method [character] A method to compute the p-value and the default value is "liu".
#'   Method "davies" represents an exact method that computes the p-value by inverting the characteristic function of the mixture chisq.
#'   Method "liu" represents an approximation method that matches the first 3 moments.
#' @param model [character]
#'
#' @return
#' * S.tau score Statistic for the variant hetergenous effect.
#' * S.pi score Statistic for the variant mean effect.
#' * p.value.S.tau P-value for testing the variant hetergenous effect.
#' * p.value.S.pi P-value for testing the variant mean effect.
#' * p.value.overall Overall p-value for testing the association between the set of SNPS/genes and outcomes.
#'   It combines p.value.S.pi and p.value.S.tau by using Fisher's procedure.
#'
#' @export
mist <- function(y, X, G, Z, method = "liu", model = c("guess", "continuous", "binary")) {
  check_y <- c("continuous", "binary")[(length(unique(y))==2)+1]
  if (any(grepl("guess", model))) {
    message('[MiST] "y" seems to be "', check_y, '", model is set to "', check_y, '"!')
    model <- check_y
  }
  if (model!=check_y) {
    warning('[MiST] "y" seems to be "', check_y,'" and model was set to "', model, '"!')
  }
  switch(
    EXPR = model,
    "continuous" = {
      message('[MiST] Linear regression is ongoing ...')
      suppressMessages(mist_linear(y = y, X = X, G = G, Z = Z, method = method))
    },
    "binary" = {
      message('[MiST] Logistic regression is ongoing ...')
      suppressMessages(mist_logit(y = y, X = X, G = G, Z = Z, method = method))
    },
    stop('[MiST] "model" must be one of "guess", "continuous" or "binary".')
  )
}


#' mist_print
#'
#' @param x [mist]
#'
#' @return list
#' @export
mist_print <- function(x) {
  cluster_name <- gsub("^GZ", "", rownames(x$out_rare))
  rownames(x$out_rare) <- NULL
  stat_rare <- cbind.data.frame(
    "SubClusters" = ifelse(any(cluster_name %in% c("", "M")), "None", cluster_name),
    x$out_rare
  )
  list(estimate = stat_rare, statistic = as.data.frame(x$out_MiST))
}


#' mist_logit
#'
#' @inheritParams mist
#'
#' @return data.frame
#' @export
mist_logit <- function(y, X, G, Z, method = "liu") {
  y <- as.vector(as.numeric(y))
  G <- as.matrix(G)
  X <- as.matrix(X)
  X <- apply(X, 2, as.numeric)
  Z <- as.matrix(Z)
  GZ <- G %*% Z
  M <- cbind(X, GZ)
  fit.0 <- stats::glm(
    formula = y ~ X - 1,
    family = stats::binomial(link = "logit")
  )
  mu.0 <- fit.0$fitted.value
  d.0 <- mu.0 * (1 - mu.0)
  res.0 <- y - mu.0
  fit.0a <- stats::glm(
    formula = y ~ -1 + X + GZ,
    family = stats::binomial(link = "logit")
  )
  mu.0a <- fit.0a$fitted.value
  d.0a <- mu.0a * (1 - mu.0a)
  res.0a <- y - mu.0a
  n <- dim(X)[1]
  I <- diag(1, n)
  D.0 <- diag(d.0)
  D.0a <- diag(d.0a)
  tXD0X <- t(X) %*% D.0 %*% X
  inv.tXD0X <- solve(tXD0X)
  tMD0aM <- t(M) %*% D.0a %*% M
  inv.tMD0aM <- solve(tMD0aM)
  P01 <- D.0 - (d.0 %o% d.0) * (X %*% (inv.tXD0X) %*% t(X))
  P02 <- D.0a - (d.0a %o% d.0a) * (M %*% (inv.tMD0aM) %*% t(M))
  S.tau <- 0.5 * t(res.0a) %*% G %*% t(G) %*% res.0a
  inv.I.pi <- solve(t(GZ) %*% P01 %*% GZ)
  S.pi <- t(res.0) %*% GZ %*% inv.I.pi %*% t(GZ) %*% res.0
  p.value.S.pi <- 1 - stats::pchisq(S.pi, df = dim(Z)[2])
  Mat <- 0.5 * t(G) %*% P02 %*% G
  eigen.value <- eigen(Mat, symmetric = TRUE)$values
  lambda <- eigen.value

  ## original code : 
  # if (method == "davies") {
  #   p.value.S.tau <- try(CompQuadForm::davies(S.tau, lambda)$Qq, silent = TRUE)
  # }
  # if (method == "liu") {
  #   p.value.S.tau <- try(c(CompQuadForm::liu(S.tau, lambda)), silent = TRUE)
  # }
  # if (class(p.value.S.tau)=="try-error") {
  #   p.value.S.tau <- NA
  # }
  ## improve tryCatch in case of error :
  if (method == "davies") {
    p.value.S.tau <- tryCatch(
      expr = CompQuadForm::davies(S.tau, lambda)$Qq,
      error = function(e) NA
    )
  }
  if (method == "liu") {
    p.value.S.tau <- tryCatch(
      expr = CompQuadForm::liu(S.tau, lambda),
      error = function(e) NA
    )
  }

  q.fisher <- -2 * (log(p.value.S.tau) + log(p.value.S.pi))
  p.value.overall <- 1 - stats::pchisq(q.fisher, df = 4)

  out_MiST <- list(
    S.pi = S.pi,
    p.value.S.pi = p.value.S.pi,
    S.tau = S.tau,
    p.value.S.tau = p.value.S.tau,
    p.value.overall = p.value.overall
  )
  
  ## extract direct effect of the cluster Pi,
  ## extract its IC, SE, P_val. 
  get_GZ <- paste0("GZ", colnames(GZ))
  try(tools::assertCondition(CI <- as.data.frame(stats::confint(fit.0a)[get_GZ, , drop = FALSE])), silent = TRUE)
  colnames(CI) <- c("CI_2.5", "CI_97.5")
  fit.0a_coef <- stats::coef(summary(fit.0a))[get_GZ, , drop = FALSE] 
  out_rare <- cbind(
    Pi_hat = fit.0a$coefficients[get_GZ],
    CI,
    SE = fit.0a_coef[, "Std. Error"],
    P_val = fit.0a_coef[, "Pr(>|z|)"], 
    OR = exp(stats::coef(fit.0a)[get_GZ])
  )
  rownames(out_rare) <- get_GZ

  list(out_MiST = out_MiST, out_rare = out_rare)
}


#' mist_linear
#'
#' @inheritParams mist
#'
#' @return data.frame
#' @export
mist_linear <- function(y, X, G, Z, method = "liu") {
  y <- as.vector(as.numeric(y))
  G <- as.matrix(G)
  X <- as.matrix(X)
  X <- apply(X, 2, as.numeric)
  Z <- as.matrix(Z)
  GZ <- G %*% Z
  M <- cbind(X, GZ)
  tXX <- t(X) %*% X
  inv.tXX <- solve(tXX)
  tMM <- t(M) %*% M
  inv.tMM <- solve(tMM)
  fit.0 <- stats::lm(y ~ X - 1)
  tilde.sigma2 <- summary(fit.0)$sigma^2
  res.0 <- fit.0$resid
  fit.0a <- stats::lm(y ~ M - 1)
  hat.sigma2 <- summary(fit.0a)$sigma^2
  res.0a <- fit.0a$resid
  n <- dim(X)[1]
  I <- diag(1, n)
  P2 <- I - X %*% inv.tXX %*% t(X)
  P1 <- I - M %*% inv.tMM %*% t(M)
  S.tau <- t(res.0a) %*% G %*% t(G) %*% res.0a
  inv.I.pi <- solve(t(GZ) %*% P2 %*% GZ)
  S.pi <- t(res.0) %*% GZ %*% inv.I.pi %*% t(GZ) %*% res.0
  S.pi <- S.pi / tilde.sigma2
  p.value.S.pi <- 1 - stats::pchisq(S.pi, df = dim(Z)[2])
  P1.G <- P1 %*% G
  Mat <- (hat.sigma2) * t(P1.G) %*% P1.G
  eigen.value <- eigen(Mat, symmetric = TRUE)$values
  lambda <- eigen.value

  ## original code : 
  # if (method == "davies") {
  #   p.value.S.tau <- try(CompQuadForm::davies(S.tau, lambda)$Qq, silent = TRUE)
  # }
  # if (method == "liu") {
  #   p.value.S.tau <- try(CompQuadForm::liu(S.tau, lambda), silent = TRUE)
  # }
  # if (class(p.value.S.tau)=="try-error") {
  #   p.value.S.tau <- NA
  # }
  ## improve tryCatch in case of error :
  if (method == "davies") {
    p.value.S.tau <- tryCatch(
      expr = CompQuadForm::davies(S.tau, lambda)$Qq,
      error = function(e) NA
    )
  }
  if (method == "liu") {
    p.value.S.tau <- tryCatch(
      expr = CompQuadForm::liu(S.tau, lambda),
      error = function(e) NA
    )
  }
  
  q.fisher <- -2 * (log(p.value.S.tau) + log(p.value.S.pi))
  p.value.overall <- 1 - stats::pchisq(q.fisher, df = 4)

  out_MiST <- list(
    S.pi = S.pi,
    p.value.S.pi = p.value.S.pi,
    S.tau = S.tau,
    p.value.S.tau = p.value.S.tau,
    p.value.overall = p.value.overall
  )

  ## extract direct effect of the cluster Pi,
  ## extract its IC, SE, P_val. 
  get_GZ <- paste0("M", colnames(GZ)) ## cf variable put in fit.0a
  try(tools::assertCondition(CI <- as.data.frame(stats::confint(fit.0a)[get_GZ, , drop = FALSE])), silent = TRUE)
  colnames(CI) <- c("CI_2.5", "CI_97.5")
  fit.0a_coef <- stats::coef(summary(fit.0a))[get_GZ, , drop = FALSE] 
  out_rare <- cbind(
    Pi_hat = fit.0a$coefficients[get_GZ],
    CI, 
    SE = fit.0a_coef[, "Std. Error"],
    P_val = fit.0a_coef[, "Pr(>|t|)"] 
  )
  rownames(out_rare) <- get_GZ

  list(out_MiST = out_MiST, out_rare = out_rare)
}
