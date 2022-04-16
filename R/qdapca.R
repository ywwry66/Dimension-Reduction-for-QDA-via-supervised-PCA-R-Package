##' This function runs the Dimension Reduction for Quadratic Discriminant Analysis via Supervised Principal Component Analysis (QDAPCA) method.
##'
##' placeholder
##' @title Dimension Reduction for Quadratic Discriminant Analysis via Supervised Principal Component Analysis
##' @param x 
##' @param y 
##' @param xnew 
##' @param rk 
##' @param include_linear 
##' @param standardize 
##' @return 
##' @author Ruiyang Wu
##' @export
qdapca <- function(x, y, xnew, rk = 1, include_linear = TRUE,
                   standardize = TRUE) {
    x0 <- x[which(y == 0), ]
    x1 <- x[which(y == 1), ]
    x0c <- scale(x0, scale = FALSE)
    x1c <- scale(x1, scale = FALSE)
    xc <- rbind(x0c, x1c)
    xc_cplx <- rbind(x0c, 1i * x1c)
    if (standardize == TRUE) {
        scaling <- diag(1 / sqrt(diag(cov(xc))))
        x <- x %*% scaling
        xnew <- xnew %*% scaling
        x0 <- x0 %*% scaling
        x1 <- x1 %*% scaling
        xc <- xc %*% scaling
        xc_cplx <- xc_cplx %*% scaling
    }
    n <- nrow(xc_cplx)
    p <- ncol(xc_cplx)
    if (n < p) {
        ei <- eigen(xc_cplx %*% t(xc_cplx))
        f <- t(xc_cplx) %*%
            ei$vectors[, sort(abs(ei$values), decreasing = TRUE,
                              index.return = TRUE)$ix[1:rk],
                       drop = FALSE]
        f <- apply(f, MARGIN = 2, FUN = to_real)
    }
    else {
        cov_diff <- t(xc_cplx) %*% xc_cplx
        mode(cov_diff) <- "double"
        ei <- eigen(cov_diff)
        f <- ei$vectors[, sort(abs(ei$values), decreasing = TRUE,
                               index.return = TRUE)$ix[1:rk],
                        drop = FALSE]
    }
    if (include_linear == TRUE) {
        m0 <- colMeans(x0)
        m1 <- colMeans(x1)
        xc_svd <- svd(xc, nu = 0)
        rank <- sum(xc_svd$d > 1e-4)
        sigma_inv <- xc_svd$v[, 1:rank] %*% diag(1 / xc_svd$d[1:rank]) %*%
            t(xc_svd$v[, 1:rank] %*% diag(1 / xc_svd$d[1:rank]))
        d <- sigma_inv %*% (m0 - m1)
        f <- cbind(f, d)
    }
    return(list(class = predict(MASS::qda(x %*% f, y), xnew %*% f)$class))
}

##' This functions fits Dimension Reduction for Quadratic Discriminant Analysis via Supervised Principal Component Analysis (QDAPCA) using K-fold cross-validation for rank/dimension
##'
##' placeholder
##' @title Fits qdapca using K-fold cross-validation for rank
##' @param x 
##' @param y 
##' @param xnew 
##' @param rk 
##' @param folds 
##' @param seed 
##' @param standardize 
##' @return 
##' @author Ruiyang Wu
##' @export
qdapca_cv <- function(x, y, xnew, rk = 1:(min(ncol(x), nrow(x), 50)),
                      folds = 5, seed = 2020,
                      include_linear = "cv",
                      standardize = TRUE) {
    if (!is.null(seed)) {
        ## reinstate system seed after simulation
        sys_seed <- .GlobalEnv$.Random.seed
        on.exit({
            if (!is.null(sys_seed)) {
                .GlobalEnv$.Random.seed <- sys_seed
            } else {
                rm(".Random.seed", envir = .GlobalEnv)
            }
        })
        set.seed(seed)
    }
    x0 <- x[which(y == 0), ]
    x1 <- x[which(y == 1), ]
    n0 <- nrow(x0)
    n1 <- nrow(x1)
    folds0 <- cut(seq_len(n0), breaks = folds, labels = FALSE)
    folds1 <- cut(seq_len(n1), breaks = folds, labels = FALSE)
    pred_err <- matrix(rep(0, 2 * length(rk)), ncol = 2)
    for (i in seq_len(folds)) {
        cat(i)
        test_ind0 <- which(folds0 == i)
        test_ind1 <- which(folds1 == i)
        test_n0 <- length(test_ind0)
        test_n1 <- length(test_ind1)
        x0_cv <- x0[-test_ind0, ]
        x1_cv <- x1[-test_ind1, ]
        x_cv <- rbind(x0_cv, x1_cv)
        xnew_cv <- rbind(x0[test_ind0, ], x1[test_ind1, ])
        y_cv <- c(rep(0, n0 - test_n0), rep(1, n1 - test_n1))
        x0c <- scale(x0_cv, scale = FALSE)
        x1c <- scale(x1_cv, scale = FALSE)
        xc <- rbind(x0c, x1c)
        xc_cplx <- rbind(x0c, 1i * x1c)
        if (standardize == TRUE) {
            scaling <- diag(1 / sqrt(diag(cov(xc))))
            x_cv <- x_cv %*% scaling
            xnew_cv <- xnew_cv %*% scaling
            x0_cv <- x0_cv %*% scaling
            x1_cv <- x1_cv %*% scaling
            xc <- xc %*% scaling
            xc_cplx <- xc_cplx %*% scaling
        }
        n <- nrow(xc_cplx)
        p <- ncol(xc_cplx)
        if (n < p)
            ei <- eigen(xc_cplx %*% t(xc_cplx))
        else {
            cov_diff <- t(xc_cplx) %*% xc_cplx
            mode(cov_diff) <- "double"
            ei <- eigen(cov_diff)
        }
        m0 <- colMeans(x0_cv)
        m1 <- colMeans(x1_cv)
        xc_svd <- svd(xc, nu = 0)
        rank <- sum(xc_svd$d > 1e-4)
        sigma_inv <- xc_svd$v[, 1:rank] %*% diag(1 / xc_svd$d[1:rank]) %*%
            t(xc_svd$v[, 1:rank] %*% diag(1 / xc_svd$d[1:rank]))
        d <- sigma_inv %*% (m0 - m1)
        for (j in seq_along(rk)) {
            if (n < p) {
                f <- t(xc_cplx) %*%
                    ei$vectors[, sort(abs(ei$values),
                                      decreasing = TRUE,
                                      index.return = TRUE)$ix[1:j],
                               drop = FALSE]
                f <- apply(f, MARGIN = 2, FUN = to_real)
            }
            else
                f <- ei$vectors[, sort(abs(ei$values),
                                       decreasing = TRUE,
                                       index.return = TRUE)$ix[1:j],
                                drop = FALSE]
            ## LDA direction not used
            if (class(try(ypred <-
                              predict(MASS::qda(x_cv %*% f, y_cv),
                                      xnew_cv %*% f)$class,
                          silent = TRUE))
                != "factor") {
                rk <- 1:(j - 1)
                break
            }
            pred_err[j, 1] <- pred_err[j, 1] +
                sum(ypred != c(rep(0, test_n0), rep(1, test_n1)))
            ## LDA direction used
            f <- cbind(f, d)
            if (class(try(ypred <-
                              predict(MASS::qda(x_cv %*% f, y_cv),
                                      xnew_cv %*% f)$class,
                          silent = TRUE))
                != "factor") {
                rk <- 1:(j - 1)
                break
            }
            pred_err[j, 2] <- pred_err[j, 2] +
                sum(ypred != c(rep(0, test_n0), rep(1, test_n1)))
        }
    }
    pred_err <- pred_err[rk, ]
    pred_err_min <- min(pred_err)
    parameter_best <- which(pred_err == pred_err_min, arr.ind = TRUE)
    rk_best <- rk[parameter_best[1, 1]]
    if (include_linear == "cv")
        include_linear <- as.logical(parameter_best[1, 2] - 1)
    out <- qdapca(x, y, xnew, rk_best, include_linear = include_linear,
                  standardize = standardize)
    return(c(out, list(rk = rk_best, include_linear = include_linear)))
}
