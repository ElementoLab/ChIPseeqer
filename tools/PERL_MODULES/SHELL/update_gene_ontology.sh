echo "getting the ontologies"
wget http://www.geneontology.org/ontology/function.ontology
wget http://www.geneontology.org/ontology/process.ontology
wget http://www.geneontology.org/ontology/component.ontology
wget http://www.geneontology.org/doc/GO.terms_and_ids
echo "getting $1 annotations"
wget http://www.geneontology.org/cgi-bin/downloadGOGA.pl/gene_association.$1.gz
gunzip gene_association.$1.gz
echo "transforming ontologies and annotations into SQL"
perl /home/olly/PERL_MODULES/SCRIPTS/gene_association2sql.pl gene_association.$1 P $1 > gene_association.$1.process.sql
perl /home/olly/PERL_MODULES/SCRIPTS/gene_association2sql.pl gene_association.$1 C $1 > gene_association.$1.component.sql
perl /home/olly/PERL_MODULES/SCRIPTS/gene_association2sql.pl gene_association.$1 F $1 > gene_association.$1.function.sql
perl ~/PERL_MODULES/SCRIPTS/terms2sql.pl < GO.terms_and_ids > GO.terms_and_ids.sql
echo "saving old data"
mysqldump $2 GO_FULL_ANNOTATION GO_CATEGORIES > dump_old.sql
echo "erasing old data"
mysql $2 -e "delete from GO_FULL_ANNOTATION; delete from GO_CATEGORIES;"
echo "loading terms"
mysql $2 -e "load data local infile 'GO.terms_and_ids.sql' into table GO_CATEGORIES;"
echo "loading data"
mysql $2 < gene_association.$1.process.sql
mysql $2 < gene_association.$1.component.sql
mysql $2 < gene_association.$1.function.sql
echo "done"
