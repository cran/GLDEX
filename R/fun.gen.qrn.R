"fun.gen.qrn" <-
function(n, dimension, scrambling,FUN="runif.sobol"){

if(FUN == "runif.sobol"){
result<-generate_sobol_set(n = n, dim = dimension,seed=scrambling)
}

if(FUN == "runif.halton"){
result<-generate_halton_faure_set(n = n+1, dim = dimension)[-1,]
}

if(FUN == "runif.sobol.owen"){
result<-generate_sobol_owen_set(n = n, dim = dimension,seed=scrambling)
}


if(FUN == "QUnif"){
result<-QUnif(n = n, p = dimension, leap = scrambling)
}

result

}

