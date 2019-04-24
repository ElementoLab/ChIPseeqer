rm /tmp/toto*
columns.pl 0 1 2 3 4 < $1 > /tmp/toto
# R script
cat <<EOF > /tmp/toto.R
a<-read.table("/tmp/toto")
a1<-c() 
a2<-c() 
for (i in 1:(dim(a)[1])) { 
 a1[i]<-sum(dbinom(a[i,4]:a[i,5], a[i,5], 0.5)); 
 a2[i]<-(1-a1[i]); 
}
a3<-cbind(a, a1, a2)
write.table(a3, file="/tmp/toto.out", sep="\t", quote=F, row.names=F, col.names=F)  
EOF
/usr/bin/R BATCH /tmp/toto.R
#cat /tmp/toto.out
perl /home/olly/PERL_MODULES/SCRIPTS/analyze_orientation_from_bestwins_output.pl /tmp/toto.out
