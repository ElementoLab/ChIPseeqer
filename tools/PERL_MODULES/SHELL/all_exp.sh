for i in `ls /home/olly/COMPARATIVE_YEAST/ORFS/BLAST_REPORTS/BINDINGMAP/SCRIPTS/OVERLAP/*.txt`;
do
echo "processing $i"
perl test_ExpressionDistribution.pl $i STRESS 
n=`basename $i`
mv test.png EXP/$n.png
done
