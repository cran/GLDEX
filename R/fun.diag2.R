"fun.diag2" <-
function(result, data, no.test = 1000,len=100,alpha=0.05)
{
data<-data[sample(length(data),no.test*len,TRUE)]
perc.rs <- rgl(len*no.test, result[1, 1], result[2, 1], result[3, 1], result[4, 1], "rs")
mm.fmkl <- rgl(len*no.test, result[1, 2], result[2, 2], result[3, 2], result[4, 2], "fmkl")
star.fmkl <- rgl(len*no.test, result[1, 3], result[2, 3], result[3, 3], result[4, 3], "fmkl")
test <- split(data, 1:no.test)
perc.rs <- split(perc.rs, 1:no.test)
mm.fmkl <- split(mm.fmkl, 1:no.test)
star.fmkl <- split(star.fmkl, 1:no.test)
rs <- sum(sapply(1:no.test, function(i,  test,  perc.rs)
ks.gof(test[[i]], perc.rs[[i]])$p.value,  test,  perc.rs) > alpha)
fmkl <- sum(sapply(1:no.test, function(i,  test, mm.fmkl)
ks.gof(test[[i]], mm.fmkl[[i]])$p.value,  test,  mm.fmkl) > alpha)
star <- sum(sapply(1:no.test, function(i,  test,  star.fmkl)
ks.gof(test[[i]], star.fmkl[[i]])$p.value,  test,  star.fmkl) > alpha)
return(cbind(rs, fmkl, star))
}

