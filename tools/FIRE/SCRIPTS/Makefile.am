bin_SCRIPTS = hypergeom

EXTRA_DIST = *.pl *.pm GroupEnrichment.pm.tmpl Hypergeom-v1.0 HEATMAPS PostScript-Simple-0.07 weblogo

if USE_MYHYPERGEOM
USE_M = 1
hypergeom : GroupEnrichment.pm
	cp GroupEnrichment.pm.tmpl GroupEnrichment.pm
	perl modify_group_enrichment_pm.pl $(USE_M)
else
USE_M = 
hypergeom : $(PWD)/Hypergeom-v1.0/hypergeom/hypergeom.c GroupEnrichment.pm
	find Hypergeom-v1.0/ -name "*.[oa]" -exec rm {} \;; cd Hypergeom-v1.0/; echo $(PWD); perl Makefile.PL PREFIX=$(PWD)/../modules; make; make install;
	cp GroupEnrichment.pm.tmpl GroupEnrichment.pm
	perl modify_group_enrichment_pm.pl $(USE_M)
endif


