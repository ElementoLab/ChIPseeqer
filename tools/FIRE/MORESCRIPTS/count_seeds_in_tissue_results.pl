BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Table;

my @T = ("721_B_lymphoblasts_721_B_lymphoblasts.txt.nodups",
"ADIPOCYTE_ADIPOCYTE.txt.nodups",
"AdrenalCortex_AdrenalCortex.txt.nodups",
"Amygdala_Amygdala.txt.nodups",
"Appendix_Appendix.txt.nodups",
"BM-CD105+Endothelial_BM-CD105+Endothelial.txt.nodups",
"BM-CD33+Myeloid_BM-CD33+Myeloid.txt.nodups",
"BM-CD34+_BM-CD34+.txt.nodups",
"BM-CD71+EarlyErythroid_BM-CD71+EarlyErythroid.txt.nodups",
"CardiacMyocytes_CardiacMyocytes.txt.nodups",
"CerebellumPeduncles_CerebellumPeduncles.txt.nodups",
"CingulateCortex_CingulateCortex.txt.nodups",
"ColorectalAdenocarcinoma_ColorectalAdenocarcinoma.txt.nodups",
"DRG_DRG.txt.nodups",
"Heart_Heart.txt.nodups",
"Hypothalamus_Hypothalamus.txt.nodups",
"Liver_Liver.txt.nodups",
"Lung_Lung.txt.nodups",
"MedullaOblongata_MedullaOblongata.txt.nodups",
"OccipitalLobe_OccipitalLobe.txt.nodups",
"OlfactoryBulb_OlfactoryBulb.txt.nodups",
"Ovary_Ovary.txt.nodups",
"PB-BDCA4+Dentritic_Cells_PB-BDCA4+Dentritic_Cells.txt.nodups",
"PB-CD14+Monocytes_PB-CD14+Monocytes.txt.nodups",
"PB-CD19+Bcells_PB-CD19+Bcells.txt.nodups",
"PB-CD4+Tcells_PB-CD4+Tcells.txt.nodups",
"PB-CD56+NKCells_PB-CD56+NKCells.txt.nodups",
"PB-CD8+Tcells_PB-CD8+Tcells.txt.nodups",
"PLACENTA_PLACENTA.txt.nodups",
"Pancreas_Pancreas.txt.nodups",
"PancreaticIslets_PancreaticIslets.txt.nodups",
"ParietalLobe_ParietalLobe.txt.nodups",
"Pituitary_Pituitary.txt.nodups",
"Pons_Pons.txt.nodups",
"PrefrontalCortex_PrefrontalCortex.txt.nodups",
"Prostate_Prostate.txt.nodups",
"SkeletalMuscle_SkeletalMuscle.txt.nodups",
"SmoothMuscle_SmoothMuscle.txt.nodups",
"SuperiorCervicalGanglion_SuperiorCervicalGanglion.txt.nodups",
"TONGUE_TONGUE.txt.nodups",
"TemporalLobe_TemporalLobe.txt.nodups",
"TestisGermCell_TestisGermCell.txt.nodups",
"TestisInterstitial_TestisInterstitial.txt.nodups",
"TestisLeydigCell_TestisLeydigCell.txt.nodups",
"TestisSeminiferousTubule_TestisSeminiferousTubule.txt.nodups",
"Thalamus_Thalamus.txt.nodups",
"Thyroid_Thyroid.txt.nodups",
"Tonsil_Tonsil.txt.nodups",
"TrigeminalGanglion_TrigeminalGanglion.txt.nodups",
"UterusCorpus_UterusCorpus.txt.nodups",
"Uterus_Uterus.txt.nodups",
"WHOLEBLOOD_WHOLEBLOOD.txt.nodups",
"WholeBrain_WholeBrain.txt.nodups",
"adrenalgland_adrenalgland.txt.nodups",
"atrioventricularnode_atrioventricularnode.txt.nodups",
"bonemarrow_bonemarrow.txt.nodups",
"bronchialepithelialcells_bronchialepithelialcells.txt.nodups",
"caudatenucleus_caudatenucleus.txt.nodups",
"cerebellum_cerebellum.txt.nodups",
"ciliaryganglion_ciliaryganglion.txt.nodups",
"fetalThyroid_fetalThyroid.txt.nodups",
"fetalbrain_fetalbrain.txt.nodups",
"fetalliver_fetalliver.txt.nodups",
"fetallung_fetallung.txt.nodups",
"globuspallidus_globuspallidus.txt.nodups",
"kidney_kidney.txt.nodups",
"leukemiachronicmyelogenousk562_leukemiachronicmyelogenousk562.txt.nodups",
"leukemialymphoblasticmolt4_leukemialymphoblasticmolt4.txt.nodups",
"leukemiapromyelocytichl60_leukemiapromyelocytichl60.txt.nodups",
"lymphnode_lymphnode.txt.nodups",
"lymphomaburkittsDaudi_lymphomaburkittsDaudi.txt.nodups",
"lymphomaburkittsRaji_lymphomaburkittsRaji.txt.nodups",
"salivarygland_salivarygland.txt.nodups",
"skin_skin.txt.nodups",
"spinalcord_spinalcord.txt.nodups",
"subthalamicnucleus_subthalamicnucleus.txt.nodups",
"testis_testis.txt.nodups",
"thymus_thymus.txt.nodups",
"trachea_trachea.txt.nodups");




my $ta = Table->new;

print "\t$ARGV[0]\n";
foreach my $r (@T) {
  
  my $f = "$r" . ".seeds7" . $ARGV[0];
  next if (! -e $f);
  $ta->loadFile($f);
  my $a_ref = $ta->getArray();
  
  my $cnt = 0;
  foreach my $r (@$a_ref) {
    if ($r->[3] >= $ARGV[1]) {
      $cnt ++;
    }
  }

  print "$r\t$cnt\n";
}
