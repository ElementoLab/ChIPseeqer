RPKM_scatterplot <- function(m, t=2) {

  xmax <- max( log(m[,1]+1) )
  ymax <- max( log(m[,2]+1) )
  
  plot(log(m[,1]+1), log(m[,2]+1), xlab=sprintf("logRPKM expression in %s", colnames(m)[1]), ylab=sprintf("logRPKM expression in %s", colnames(m)[2]), xlim=c(0,xmax*1.05), ylim=c(0,ymax*1.05))
  
  ra <- log(m[,1]+1) - log(m[,2]+1) 
  do <- ra > log2(t)
  up <- ra < log2(1/t)
  par(new=T)
  plot(log(m[do,1]+1), log(m[do,2]+1), xlab=sprintf("logRPKM expression in %s", colnames(m)[1]), ylab=sprintf("logRPKM expression in %s", colnames(m)[2]), xlim=c(0,xmax*1.05), ylim=c(0,ymax*1.05), col="red")

    par(new=T)
   plot(log(m[up,1]+1), log(m[up,2]+1), xlab=sprintf("logRPKM expression in %s", colnames(m)[1]), ylab=sprintf("logRPKM expression in %s", colnames(m)[2]), xlim=c(0,xmax*1.05), ylim=c(0,ymax*1.05), col="green")

  abline(log2(t), 1)
  abline(log2(1/t), 1)

}