#mCP  <- read.table("bcl6peaks.txt.RefGene.DISTPEAKS.2kbaround.cons", row.names=1)
#mRP  <- read.table("bcl6peaks.txt.RefGene.DISTPEAKS.2kbaround.randcons", row.names=1)

#mCP  <- read.table("toto", row.names=1)
#mRP  <- read.table("tata", row.names=1)

plotAverageReadProfile <- function(fCP, fRP=NA) {

mCP  <- read.table(fCP, row.names=1)

nCP  <- dim(mCP)[1]
nb   <- dim(mCP)[2]
ymin <-  min(colMeans(mCP))
ymax <-  max(colMeans(mCP))
par(cex=1.3)
plot(1:nb, colMeans(mCP), type="l", lwd=5, col="red", ylab="Average Read Density", xlab="Distance to motif (bp)", ylim=c(ymin, ymax), main=fCP) 

  
if (!is.na(fRP)) {
 mRP  <- read.table(fRP, row.names=1)
 par(new=T)
 plot(10*(1:nb)-1000, colMeans(mRP), type="l", lwd=5, col="green", ylab="", xlab="", ylim=c(ymin, ymax))
 cleg <- sprintf("%d distal peaks", nCP)
 rleg <- sprintf("%d random peaks", nRP)
 legend("topleft", c(cleg, rleg), col=c("red", "green"), lty=1, lwd=5)
}

}



if (! interactive()){
        # invoked from command line
        i= commandArgs(TRUE)
        if (length(i)!= 1){
                writeLines("Error!\nUse R --slave --args FILE < makeAvgConsProfile.R \n") 
                quit(save='no')
        }

        pdffile <- sprintf("%s.pdf", i[1])
 	pdf(file=pdffile)
        plotAverageReadProfile(i[1]);        


	dev.off()
	system(sprintf("uuencode %s %s | mail -s \"%s\" ole2001@med.cornell.edu ", pdffile, pdffile, i[1]))
}
