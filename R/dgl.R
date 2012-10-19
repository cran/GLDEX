"dgl" <-
function (x, lambda1 = 0, lambda2 = NULL, lambda3 = NULL, lambda4 = NULL, 
    param = "fmkl", inverse.eps =1e-08, max.iterations = 500) 
{
    lambdas <- .gl.parameter.tidy(lambda1, lambda2, lambda3, 
        lambda4, param)
    if (!gl.check.lambda.alt1(lambdas, param = param, vect = TRUE)) {
         return(rep(NA,length(x)))
    }
    extreme <- qgl(c(0, 1), lambda1 = lambdas, param = param)
    outside.range <- !as.logical((x <= extreme[2]) * (x >= extreme[1]))
    u <- pgl(x, lambdas, param = param, inverse.eps = inverse.eps, 
        max.iterations = max.iterations)
    dens <- qdgl(u, lambda1 = lambdas, param = param)
    dens[outside.range] <- 0
    dens
}

