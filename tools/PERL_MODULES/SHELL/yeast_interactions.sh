#!/bin/bash
o=`yeast_orfname.pl $1`
echo $o
if [[ $o == "" ]]
then 
o=$1
fi
mysql YEASTHOUSE -e "select INTERACTIONS.*, YEASTHOUSE.GENENAME, GO_ANNOTATION.DESCRIPTION from INTERACTIONS,YEASTHOUSE, GO_ANNOTATION where P1 = '$o' and YEASTHOUSE.ORF = INTERACTIONS.P2 and GO_ANNOTATION.ORF = INTERACTIONS.P2;"

