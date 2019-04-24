for f in `ls /home/olly/DATA/YEASTS/EXPRESSION/CLEANED/*.clean`;
do
n=`basename $f`
/home/olly/PROGRAMS/COEXPRESSION/calculate_expression_coherence_for_all_sets -expfile $f -listfile $1 > $n.correlations
done
