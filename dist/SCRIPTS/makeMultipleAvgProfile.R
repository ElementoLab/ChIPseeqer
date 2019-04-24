
plotMultipleAverageProfile <- function(a_f) {

 # load data 
 m    <- c()
 ymin <- 100
 ymax <- -100
 for (f in a_f) {
   mCP  <- read.table(f, row.names=1, header=F)
   ymin <- min(ymin, colMeans(mCP))
   ymax <- max(ymax, colMeans(mCP))
   m    <- rbind(m, colMeans(mCP))
 }

 par(cex=1.3)
 n  <- dim(m)[1]
 nb <- dim(m)[2]
 col <- c("red", "green", "blue", "black", "orange", "lightblue", "yellow", "purple")
 for (i in 1:n) {
   plot(10*(1:nb)-1000, m[i,], type="l", lwd=5, col=col[i], ylab="Reads", xlab="Distance to TSS (bp)", ylim=c(ymin, ymax)) 
   par(new=T)
 }
 
 a_f <- sub("^.*profiles.txt", "", a_f)
 legend("topleft", basename(a_f), col=col[1:n], lty=1, lwd=5)

}


if (! interactive()){
        # invoked from command line
        i= commandArgs(TRUE)
        if (length(i) == 1){
                writeLines("Error!\nUse R --slave --args FILE < makeAvgConsProfile.R \n") 
                quit(save='no')
        }

        pdffile <- "out.pdf"
 	pdf(file=pdffile)
        plotMultipleAverageProfile(i)        


	dev.off()
	system(sprintf("uuencode %s %s | mail ole2001@med.cornell.edu ", pdffile, pdffile))
}
