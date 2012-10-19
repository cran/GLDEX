"pgl" <-
function (q, lambda1 = 0, lambda2 = NULL, lambda3 = NULL, lambda4 = NULL, 
    param = "fmkl", inverse.eps = 1e-08, max.iterations = 500) 
{
    lambdas <- .gl.parameter.tidy(lambda1, lambda2, lambda3, 
        lambda4, param)
    if (!gl.check.lambda.alt1(lambdas, param = param, vect = TRUE)) {
       return(rep(NA,length(q)))
    }
    jr <- q
    jr[sort.list(q)] <- seq(along = q)
    order.x <- order(q)
    xx <- sort(q)
    extreme <- qgl(c(inverse.eps, 1 - inverse.eps), lambda1 = lambdas, 
        param = param)
    max.e <- extreme[2]
    min.e <- extreme[1]
    ind.min <- xx <= min.e
    ind.max <- xx >= max.e
    q <- xx[as.logical((xx < max.e) * (xx > min.e))]
    length.of.vector <- length(q)
    u <- 0 * q
    result <- switch(param, freimer = , frm = , FMKL = , fmkl = .C("gl_fmkl_distfunc", 
        lambdas[1], lambdas[2], lambdas[3], lambdas[4], as.double(0), 
        as.double(1), inverse.eps, as.integer(max.iterations), 
        as.double(q), as.double(u), as.integer(length.of.vector), 
        PACKAGE = "GLDEX"), ramberg = , ram = , RS = , rs = .C("gl_rs_distfunc", 
        lambdas[1], lambdas[2], lambdas[3], lambdas[4], as.double(0), 
        as.double(1), inverse.eps, as.integer(max.iterations), 
        as.double(q), as.double(u), as.integer(length.of.vector), 
        PACKAGE = "GLDEX"), stop("Error: Parameterisation must be fmkl or rs"))
    if (!(is.numeric(result[[1]]))) {
        stop("Values for quantiles outside range. This shouldn't happen")
    }

    else {
        u <- result[[10]]
    }
    xx[as.logical((xx < max.e) * (xx > min.e))] <- u
    xx[ind.min] <- 0
    xx[ind.max] <- 1
    xx[jr]
}

