
plotMultipleAverageProfile <- function(a_f, genelistfile=NA, title=NA, ws=10) {

 # load genelist
 if (!is.na(genelistfile)) {
   d    <- as.vector(read.table(genelistfile, header=F)[,1])
 }

 # load data
 m    <<-  c()
 ymin <-  100
 ymax <- -100
 for (f in a_f) {
   print(sprintf("Loading %s", f))
   mCP  <- read.table(f, row.names=1, header=F)
   print("Done")
   rome <- rowMeans(mCP)  # mean for all rows
   med  <- mean(rome) 
   print(med)
   mCP  <- mCP * 5 / med  # set median target to 5
   if (!is.na(genelistfile)) { 
     dCP  <- mCP[d,]        # get
     dCP  <- dCP[!is.na(dCP[,1]),] # should get rid of all NA 
     print(dim(dCP))
   } else {
     dCP <- mCP
   }
   ymin <- min(ymin, colMeans(dCP))
   ymax <- max(ymax, colMeans(dCP))
   m    <<- rbind(m, colMeans(dCP))
 }


 par(cex=1.3)
 n  <- dim(m)[1]
 nb <- dim(m)[2]
  print(nb)
 col <- c("green", "red", "blue", "black", "orange", "lightblue", "yellow", "purple")
 for (i in 1:n) {
   plot(ws*(1:nb)-nb*ws/2, m[i,], type="l", lwd=5, col=col[i], ylab="Reads", xlab="Distance to TSS (bp)", ylim=c(ymin, ymax), main=title) 
   par(new=T)
 }
 par(new=F) 
 a_f <- sub("^.+(NUCL.+)$", "\\1", a_f)
# legend("topleft", basename(dirname(a_f)), col=col[1:n], lty=1, lwd=5)
  legend("topleft", a_f, col=col[1:n], lty=1, lwd=5)

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
 
        genelist <- i[1]
        i        <- i[-1]
        plotMultipleAverageProfile(i, genelist)        


	dev.off()
	system(sprintf("uuencode %s %s | mail ole2001@med.cornell.edu ", pdffile, pdffile))
}
