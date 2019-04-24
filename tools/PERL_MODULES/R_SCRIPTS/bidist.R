bidist<-function(a, b, bks, ymax, leg=c("actual", "randomized"), xleg=15) {
ab<-c(a,b)
hab<-hist(ab, br=bks)
br<-hab$breaks
mi<-hab$mids

ha<-hist(a, br=br)
da<-ha$density
da<-da/sum(da)

db<-hist(b, br=br)$density
db<-db/sum(db)
dd<-data.frame(da, db)

mix<-c()
nb<-dim(t(mi))[2]

for (j in 1:nb) {
 mix[j]<-NA
}
for (i in 0:10) {
 n<-i*50;
 best_j<- 0
 best_s<- 500000
 for (j in 1:nb) {
  if (((mi[j]-n)>=0) & ((mi[j]-n)<best_s)) {
    best_s<-mi[j]-n
    best_j<-j
  }
 }
 if (best_j <= nb) {
   mix[best_j]<-n
 }
}
print(mix);

par(cex=1.5)
#barplot(t(dd), beside=T,  ylim=c(0.0, max(max(da),max(db))), col=c("red", "black"), names.arg=mi, border=F)

barplot(t(dd), beside=T,  ylim=c(0.0, ymax), col=c("red", "black"), names.arg=mi, border=F)

par(cex=1.2)
#legend(0, max(max(da),max(db)), leg, c("red", "black"))

legend(xleg, ymax, leg, c("red", "black"))

#barplot(t(dd), beside=T,  ylim=c(0.0, ymax), col=c("red", "black"), names.arg=mix, border=F, xpd=F, legend.text=c("actual", "randomized"))
#barplot(t(dd), beside=T,  ylim=c(0.0, max(max(da),max(db))), col=c("red", "black"), names.arg=mix, border=F, xpd=F)

mi
}


tridist<-function(a, b, c, bks) {
abc<-c(a,b, c)
habc<-hist(abc, br=bks)
br<-habc$breaks
mi<-habc$mids

ha<-hist(a, br=br)
da<-ha$density
da<-da/sum(da)

db<-hist(b, br=br)$density
db<-db/sum(db)

dc<-hist(c, br=br)$density
dc<-dc/sum(dc)


dd<-data.frame(da, db, dc)



par(cex=1.5)
barplot(t(dd), beside=T,  ylim=c(0.0, 0.25), col=c("red", "black", "yellow"), names.arg=mi, border=F, xpd=F, legend.text=c("conserved", "non-conserved", "random"))
#barplot(t(dd), beside=T,  ylim=c(0.0, max(max(da),max(db), max(dc))), col=c("red", "black", "yellow"), names.arg=mi, border=F, xpd=F, legend.text=c("conserved", "non-conserved", "random"))

#mi
}
