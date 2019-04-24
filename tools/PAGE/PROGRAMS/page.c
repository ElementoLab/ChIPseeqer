#define _GNU_SOURCE
#include <search.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <search.h>
#include <limits.h>
#include <time.h>
#include "statistics.h"

#ifdef USEMPATROL
#include <mpatrol.h>
#endif


#include "dataio.h"
#include "statistics.h"
#include "information.h"
#include "mi_library.h"
#include "hashtable.h"



#define min(a, b) ((a) < (b) ? a : b)
#define max(a, b) ((a) > (b) ? a : b)

#define MIN_KMER_COUNT 20

int CmpFunc(const void* _a, const void* _b) ;
void read_go_index_data(char *filename, int *gene_num, int **line_sizes, char**** go_terms, char*** gene_names, char*** go_cat_list, int go_terms_num, int cat_min_count, int cat_max_count, struct my_hsearch_data *hash_excluded) ;
void read_go_names_data(char *filename, int *go_terms_num, char*** go_terms, char*** go_desc, char* cat_status, struct my_hsearch_data *hash_excluded) ;
double minCondInfoNormalized(int *A_go_indices, unsigned char** go_profile, int DA, int nA, int* B_q, int DB, int* E_q, int DE, int n, int shuffle, double minr, int v, int* midx) ;
double evalSeed(int ii, unsigned char** gene_kmer, int nborfs, double score, int mbins, int* E_q, int ebins, int shuffle, int jn, int jn_f, int jn_t) ;
int get_cat_id(char c) ;
double my_max_rank_test(double score, int* M_q, int mbins, int* E_q, int ebins, int nborfs, int shuffle) ;

char**   ch;
char**   d_ch;
int      nbchars = 15;
float**  p_ch;



typedef struct _GO {
  int    index;      // where is it found
  double score;
  double pearson;
} GO;



int main(int argc, char ** argv)
{

 // GENERAL
 int i, j, nk;

 int divbins = 50 ;
 int *A_go_indices ;
 int *B_go_indices ;
 int nb_cat = 0;
 int gene_num ;
 int *line_sizes ;
 char ***go_terms ;
 char **gene_names ;
 struct my_hsearch_data *hash_go_index ;
 int go_terms_num ;
 char **go_terms_list ;
 char **go_desc ;
 struct my_hsearch_data *hash_go_names ;
 struct my_hsearch_data *hash_go_excluded ;
 unsigned char **go_profile ;
 GO *go_array ;
 int cluster_num = 0 ;
 int *cluster_gene_count ;
 char cat_status[3] ;
 int cat_max_count = 100000 ;
 int cat_min_count = 0  ;
 int independence=1 ;

 float *pvalues ;
 float max_p = 0.005;

 int* line_sizes_filtered;
 char*** go_terms_filtered;

 double **pvalue ;
 double** pvalue_u;

 ENTRY e;
 ENTRY *ep;

 char*   goindexfile ;
 char*   gonamesfile ;
 char*   pvaluematrixfile ;
 char*   countmatrixfile;

 FILE*   f ;
 FILE*   fc;

 char*   expfile;
 int     m;
 char**  rownames;
 char**  tmprownames;
 char**  colnames;
 float** data;

 int max_seq = 1000000 ;
 int*     idx_gene_name;
 float*   E;
 int*     E_q;
 float*   E_q_bins = 0;
 int*     M_q;

 int      quantized = 1;

 int      ebins = 0;
 int      mbins   = 2;

 int      shuffle = 10000;

 float    minr    = 5.0; // min imi for accepting a motif
 double   minratio;

 int      midx;
 char*    logfile = 0;
 FILE*    fplog = 0;
 int      jn   = 0;
 int      jn_t = 0;
 int      jn_f = 3;
 int      idxgo;
 int      idx = 0;
 char*    go;
 int      idj;

 int      hashret = 0;
 FILE     *kill ;
 char     killprocess[500] ;
 char*    pathwaylistfile = 0;
 int      n_pathwaylist = 0;
 char**   pathwaylist = 0;
 struct my_hsearch_data *hash_pathwaylist;
 
 int      correl = 0;
 
 char*    outprofiles  = 0;
 char*    outquantized = 0;
 
 int      verbose = 0;

 if (argc == 1) {
   printf("Usage : page -goindexfile FILE -gonamesfile FILE -expfile FILE -quantized 0/1\n");
   exit(0);
 }

 cat_status[get_cat_id('F')] = 1 ;
 cat_status[get_cat_id('P')] = 1 ;
 cat_status[get_cat_id('C')] = 1 ;
 
 expfile                =      get_parameter(argc, argv, "-expfile");
 goindexfile            =      get_parameter(argc, argv, "-goindexfile");
 gonamesfile            =      get_parameter(argc, argv, "-gonamesfile");
 quantized              = atoi(get_parameter(argc, argv, "-quantized"));
 pvaluematrixfile       =      get_parameter(argc, argv, "-pvaluematrixfile");
 countmatrixfile        =      get_parameter(argc, argv, "-countmatrixfile");
 if (exist_parameter(argc, argv, "-logfile")) {
   logfile         = get_parameter(argc, argv, "-logfile");
   
   fplog = fopen(logfile, "w");
   if (!fplog) {
     die("Cannot open logfile\n");
   }
   
 }
 
 if (exist_parameter(argc, argv, "-minr"))
   minr = atof(get_parameter(argc, argv, "-minr")) ;
 
 if (exist_parameter(argc, argv, "-jn_t"))
   jn_t = atoi(get_parameter(argc, argv, "-jn_t")) ;

 if (exist_parameter(argc, argv, "-correl"))
   correl = atoi(get_parameter(argc, argv, "-correl")) ;

 if (exist_parameter(argc, argv, "-verbose"))
   verbose = atoi(get_parameter(argc, argv, "-verbose")) ;

 if (exist_parameter(argc, argv, "-ebins"))
   ebins = atoi(get_parameter(argc, argv, "-ebins")) ;

 if (exist_parameter(argc, argv, "-divbins"))
   divbins = atoi(get_parameter(argc, argv, "-divbins")) ;
 
 if (exist_parameter(argc, argv, "-shuffle"))
   shuffle = atoi(get_parameter(argc, argv, "-shuffle")) ;
 
 if (exist_parameter(argc, argv, "-F"))
   cat_status[get_cat_id('F')] = atoi(get_parameter(argc, argv, "-F")) ;

 if (exist_parameter(argc, argv, "-P"))
   cat_status[get_cat_id('P')] = atoi(get_parameter(argc, argv, "-P")) ;

 if (exist_parameter(argc, argv, "-C"))
   cat_status[get_cat_id('C')] = atoi(get_parameter(argc, argv, "-C")) ;

 if (exist_parameter(argc, argv, "-catmaxcount"))
   cat_max_count = atoi(get_parameter(argc, argv, "-catmaxcount")) ;

 if (exist_parameter(argc, argv, "-catmincount"))
   cat_min_count = atoi(get_parameter(argc, argv, "-catmincount")) ;

 if (exist_parameter(argc, argv, "-independence"))
   independence = atoi(get_parameter(argc, argv, "-independence")) ;

 if (exist_parameter(argc, argv, "-max_p"))
   max_p = atof(get_parameter(argc, argv, "-max_p")) ;

 if (exist_parameter(argc, argv, "-pathwaylist"))
   pathwaylistfile = get_parameter(argc, argv, "-pathwaylist") ;

 if (exist_parameter(argc, argv, "-outprofiles"))
   outprofiles = get_parameter(argc, argv, "-outprofiles") ;

 if (exist_parameter(argc, argv, "-outquantized"))
   outquantized = get_parameter(argc, argv, "-outquantized") ;
 

 strcpy (killprocess, expfile) ;
 strcat (killprocess, ".kill") ;

 //
 //  read in expression data
 //
 printf("Reading expression data ... ");
 readFloatTable(expfile, &m, &max_seq, &data, &rownames, &colnames, 0, 1);
 printf("Done.\n");

 fflush(stdout);

 hash_go_excluded = (struct my_hsearch_data*)calloc(1,  sizeof(struct my_hsearch_data));

 printf("Creating hash ... "); fflush(stdout);
 hashret = my_hcreate_r(100000, hash_go_excluded);
 if (hashret == 0) {
   printf("Could not create hash table ...\n");
   exit(0);
 }
 printf("Done\n");
 kill = NULL ;
 kill = fopen (killprocess, "r") ;
 if(kill != NULL){
   exit(0) ;
 }


 //
 // read GO names table
 //
 printf("Reading pathway names table ... ");

 // hash_go_excluded contains a list of all pathways that we excluded, for diverse reasons (P,C,F .. or not in pathway list)
 read_go_names_data(gonamesfile, &go_terms_num, &go_terms_list, &go_desc, cat_status, hash_go_excluded) ;
 printf("Done.\n");
 kill = NULL ;
 kill = fopen (killprocess, "r") ;
 if(kill != NULL){
   exit(0) ;
 }


 //
 // in case there is a list of pathways to keep
 //
 if (pathwaylistfile != 0) {
   
   pathwaylist = readStringSet(pathwaylistfile, &n_pathwaylist); 
   
   // build a hash table 
   hash_pathwaylist = (struct my_hsearch_data*)calloc(1,  sizeof(struct my_hsearch_data));

   hashret = my_hcreate_r(100000, hash_pathwaylist);
   if (hashret == 0) {
     printf("Could not create hash table ...\n");
     exit(0);
   }
   
   if (verbose == 1) {
     printf("Considering only %d pathways.\n", n_pathwaylist);
   }


   // add 
   for  (i=0 ; i<n_pathwaylist ; i++) {
     e.key   = strdup(pathwaylist[i]);
     e.data  = (char*)i;
     hashret = my_hsearch_r(e, ENTER, &ep, hash_pathwaylist);
     if (hashret == 0) {
       printf("Could not enter entry into hash table ...\n");
       exit(0);
     }     
   } // end for


   // now traverse the whole list of pathways to exclude the above
   for (i=0 ; i<go_terms_num ; i++) {
     
     // look it up in pathwaylist   
     e.key = go_terms_list[i];        
     my_hsearch_r(e, FIND, &ep, hash_pathwaylist) ;
     if (!ep) {
       
       //
       // if not there, add it to exclusion list      
       //
       e.key = strdup(go_terms_list[i]) ;
       e.data = (char*) 1 ;
       hashret = my_hsearch_r(e, ENTER, &ep, hash_go_excluded);
       if (hashret == 0) {
	 printf("Could not enter entry into hash_excluded ...\n");
	 exit(0);
       }
       
     }

 
   }


 } // end if


 //
 // read GO index table
 //
 printf("Reading pathway index table ... ");

 // 
 // returns 
 //     go_terms   (list of go terms for gene)
 //     line_sizes (# GO terms)
 //      
 // 
 read_go_index_data(goindexfile, &gene_num, &line_sizes, &go_terms, &gene_names, 
		    &go_terms_list, go_terms_num, cat_min_count, cat_max_count, hash_go_excluded) ;

 printf("Done.\n");

 printf("Number of pathways read: %d\n", go_terms_num) ;

 fflush(stdout);
 kill = NULL ;
 kill = fopen (killprocess, "r") ;
 if(kill != NULL){
   exit(0) ;
 }

 //
 //  create a hash for gene names in gene ontology index file
 //
 hash_go_index = (struct my_hsearch_data*)calloc(1,  sizeof(struct my_hsearch_data));

 printf("Creating hash table for gene names in pathway index file ... ");
 hashret = my_hcreate_r(100000, hash_go_index);
 if (hashret == 0) {
   printf("Could not create hash table ...\n");
   exit(0);
 }
 for (i=0 ; i<gene_num ; i++)
 {
   e.key = strdup(gene_names[i]);
   e.data = (char*)i;
   hashret = my_hsearch_r(e, ENTER, &ep, hash_go_index);
   if (hashret == 0) {
     printf("Could not enter entry into hash table (table full?) ...\n");
     exit(0);
   }
 }
 printf("Done.\n");
 kill = NULL ;
 kill = fopen (killprocess, "r") ;
 if(kill != NULL){
   exit(0) ;
 }


 idx_gene_name = (int*)  malloc( max_seq * sizeof(int));
 E             = (float*)malloc( max_seq * sizeof(int));



 //
 //   make sure we retain only genes BOTH in expression file AND go index file
 //

 line_sizes_filtered = (int*)   malloc(max_seq * sizeof(int));
 go_terms_filtered   = (char***)malloc(max_seq * sizeof(char**));
 tmprownames         =  (char**)malloc(max_seq * sizeof(char*));

 idx = 0;
 for (i=0; i<max_seq; i++) {

   e.key = rownames[i];
   my_hsearch_r(e, FIND, &ep, hash_go_index);

   if (ep) {
     idxgo                    = (int)(ep->data);
     E [idx]                  = data[i][0];
     tmprownames[idx]         = rownames[i];
     line_sizes_filtered[idx] = line_sizes[idxgo];
     go_terms_filtered[idx]   = go_terms[idxgo];
     idx ++;
   }
 
 }

 rownames           = tmprownames;
 line_sizes         = line_sizes_filtered;
 go_terms           = go_terms_filtered;
 gene_num           = idx;

 //
 //  create a hash index for gene ontology names; "pathway ID" => idx
 //

 hash_go_names = (struct my_hsearch_data*)calloc(1,  sizeof(struct my_hsearch_data));
 my_hcreate_r(100000, hash_go_names);

 // go_terms_num was read in read_go_names_data
 for (i=0 ; i<go_terms_num ; i++) {
   e.key = strdup(go_terms_list[i]);
   e.data = (char*)i;
   hashret = my_hsearch_r(e, ENTER, &ep, hash_go_names);
   if (hashret == 0) {
     printf("Could not enter entry into hash table ...\n");
     exit(0);
   }
 }
 printf("Done.\n");


 //
 // Populate the go profile table
 //
 printf("Creating a table of GO profile, # genes = %d by # terms = %d ... \n", gene_num, go_terms_num) ;

 kill = NULL ;
 kill = fopen (killprocess, "r") ;
 if(kill != NULL){
   exit(0) ;
 }

 go_profile = (unsigned char**)malloc (gene_num * sizeof (unsigned char*)) ;
 for (i=0 ; i<gene_num ; i++) {
   go_profile[i] = create_binarized_array(go_terms_num);
 }


 //
 // there is something to change here ... 
 // because pathways that are excluded still have a col
 // full of zeros, and their MI is evaluated 
 //
 for (i=0 ; i<gene_num ; i++) {

   for (j=0 ; j<line_sizes[i] ; j++) {

     go    = go_terms[i][j] ;
     e.key = go;
     my_hsearch_r(e, FIND, &ep, hash_go_names) ;
     if (!ep) {
       continue ;
     }

     idj = (int)ep->data ;

     set_entry_in_binarized_array(go_profile[i], idj);

   }
 }
 printf("Done.\n");
 kill = NULL ;
 kill = fopen (killprocess, "r") ;
 if(kill != NULL){
   exit(0) ;
 }

 my_hdestroy_r(hash_go_names);
 //my_hdestroy_r(hash_go_excluded);

 //
 //  determine number of bins
 //
 if ((quantized == 0) && (ebins == 0)) {
   ebins = (int)(0.5 + (float)gene_num / ( divbins * mbins ));
 }


 if (quantized == 0) {
   // add a little random number to each value in the E vector
   add_small_values_to_identical_floats(E, gene_num);
 }

 //
 //  quantize expression profile
 //
 quantize_E(E, gene_num, quantized, &ebins, &E_q, &E_q_bins);


 if ((quantized == 0) && (outquantized != 0)) {
   FILE* fpoq = fopen(outquantized, "w");
   if (!fpoq) {
     printf("Cannot open outquantized file %s\n", outquantized);
     exit(0);
   }
   fprintf(fpoq, "GENE\tINDEX\n");
   for (i=0; i<gene_num; i++) {
     fprintf(fpoq, "%s\t%d\n", rownames[i], E_q[i]);
   }
   fclose(fpoq);
 }


 go_array = (GO*)malloc(go_terms_num * sizeof(GO)) ;

 printf("go_terms_num = %d\n", go_terms_num);


 for (nk=0; nk<go_terms_num; nk++) {
   kill = NULL ;
   kill = fopen (killprocess, "r") ;
   if(kill != NULL){
     exit(0) ;
   }


   // COMPUTE initial profile
   M_q = c_matrix_column_binarized(go_profile, nk, gene_num);

   //  calculate mi value
   float mi = CalculateMIbasic(M_q, E_q, gene_num, mbins, ebins);

   // store
   go_array[nk].index = nk;
   go_array[nk].score = mi;

   if (quantized == 0) {

     // calculate correlation (in case 
     float* M_q_continuous = itof_vector(M_q, gene_num);
     float  pe = pearson(M_q_continuous, E, gene_num);     
     go_array[nk].pearson = pe;

     free(M_q_continuous);

   }

   free(M_q);

 }

 qsort((void*)go_array, go_terms_num, sizeof(GO), CmpFunc);

 for (nk=0; nk<min(15,go_terms_num); nk++) {
   printf("cat %d, mi = %4.3f, %s\n", go_array[nk].index, go_array[nk].score, go_desc[ go_array[nk].index ] );
   fflush(stdout);
 }
 fflush(stdout);
 printf("\nRESULTS\n");

 B_go_indices  = (int*)malloc(go_terms_num * sizeof(int));
 A_go_indices  = (int*)malloc(go_terms_num * sizeof(int));

 nb_cat = 0;

 int cnt_notpassed = 0;
 pvalues = (float*)malloc(go_terms_num * sizeof(float));

 //
 // start the testing loop
 //
 for (nk=0; nk<go_terms_num; nk++) {
   

   kill = NULL ;
   kill = fopen (killprocess, "r") ;
   if(kill != NULL){
     exit(0) ;
   }

   // look up whether this pathway is in the excluded 
   e.key = go_terms_list[ go_array[nk].index ];        
   my_hsearch_r(e, FIND, &ep, hash_go_excluded);

   if (ep) {
     continue;
   }
 
   // if we only want to show positive correl
   if ((correl == 1) && (go_array[nk].pearson < -1e-10)) {
     continue;
   }
   
   // if we only want to show negative correl
   if ((correl == -1) && (go_array[nk].pearson > 1e-10)) {
     continue;
   }


   //printf("Evaluating cat %d, mi = %4.3f, %s\n", go_array[nk].index, go_array[nk].score, go_desc[ go_array[nk].index ] );

   if (cnt_notpassed == 20) {
     printf("Total number of categories evaluated = %d.\n", nk - 10);
     fflush(stdout);
     break;
   }

   //
   //  conditional information high enough ?
   //

   M_q = c_matrix_column_binarized(go_profile, go_array[nk].index, gene_num);

   if (nb_cat > 0) {
     minratio = minCondInfoNormalized(A_go_indices, go_profile, mbins, nb_cat, M_q, mbins, E_q, ebins, gene_num, shuffle, minr, 1, &midx);
   } else {
     minratio = minr+1;
   }

   free(M_q);





   if ((minratio < minr) && (independence == 1))
     {
     //printf("Did not minr test (%4.3f)\n", minratio);
     //printf("not optimized, minratio = %4.3f (< minr=%4.3f).\n", minratio, minr);

     if (logfile != 0) {

       //if (pass > max_p) {  // login only if passes
	 fprintf(fplog, "[%s,%s] was killed by [%s,%s]\n", go_terms_list[ go_array[nk].index ], go_desc[ go_array[nk].index ], go_terms_list[ A_go_indices[midx] ], go_desc[ A_go_indices[midx] ]);
	 //}

     }

     continue;

   }

   //
   //  statistical test (shufflings + jack-knife)
   //
   double pass = evalSeed(go_array[nk].index, go_profile, gene_num, go_array[nk].score, mbins, E_q, ebins, shuffle, jn, jn_f, jn_t);


   if (pass >= max(0.999,max_p)) {
     //printf("Did not pass stat test\n");
     cnt_notpassed ++;
     printf("Didn't pass: cat %d, mi = %4.3f, pvalue = %5.4f, %s\n", go_array[nk].index, go_array[nk].score, pass, go_desc[ go_array[nk].index ] );
     fflush(stdout);
     continue;
   }

   if (pass >max_p){
     cnt_notpassed ++;
     printf("Didn't pass (cutoff=%1.4f): cat %d, mi = %4.3f, pvalue = %5.4f, %s\n", max_p, go_array[nk].index, go_array[nk].score, pass, go_desc[ go_array[nk].index ] );
     fflush(stdout);
     continue;
   }




   cnt_notpassed = 0;



   A_go_indices[ nb_cat ] = go_array[nk].index;
   pvalues[ nb_cat ] = pass ;
   nb_cat++ ;

   printf("Passed: cat %d, mi = %4.3f, pvalue = %5.4f, %s\n", go_array[nk].index, go_array[nk].score, pass, go_desc[ go_array[nk].index ] );
   fflush(stdout);
 }

 printf("Number of categories that passed the tests = %d\n", nb_cat) ;
 fflush(stdout);

 //
 // defining the number of clusters
 //
 for (i=0 ; i<gene_num ; i++){
   if (E_q[i]>cluster_num){
     cluster_num = E_q[i] ;
   }
 }
 cluster_num++ ;
 printf("Number of clusters: %d\n", cluster_num) ;

 //
 // counting the number of genes in each cluster
 //
 cluster_gene_count = (int*)malloc(cluster_num * sizeof(int)) ;
 for (i=0 ; i<cluster_num ; i++){
   cluster_gene_count[i] = 0 ;
 }
 for (i=0 ; i<gene_num ; i++){
   cluster_gene_count[E_q[i]]++ ;
 }

 //
 //calculating p-values
 //
 printf("calculating the p-value matrix.\n") ;


 printf("Allocating memory ... ");
 pvalue = (double**)malloc(cluster_num*sizeof(double*)) ;
 for (i=0 ; i<cluster_num ; i++){
   pvalue[i] = (double*)malloc(nb_cat*sizeof(double)) ;
 }

 pvalue_u = (double**)malloc(cluster_num*sizeof(double*)) ;
 for (i=0 ; i<cluster_num ; i++){
   pvalue_u[i] = (double*)malloc(nb_cat*sizeof(double)) ;
 }
 printf("Done\n");

 int** gocat_cluster_count;
 gocat_cluster_count = (int**)malloc( nb_cat * sizeof(int*));
 for (i=0; i<nb_cat; i++) {
   gocat_cluster_count[i] = (int*)calloc( cluster_num, sizeof(int));
 }


 for (i=0 ; i<nb_cat ; i++){
   int go_index = A_go_indices[i] ;
   M_q = c_matrix_column_binarized(go_profile, go_index, gene_num);

   int *cluster_go_count ; //count the go cat in each cluster
   int cat_gene_count=0 ; //count the number of genes in the go cat
   cluster_go_count = (int*)malloc(cluster_num * sizeof(int)) ;
   for (j=0 ; j<cluster_num ; j++){
     cluster_go_count[j]=0 ;
     gocat_cluster_count[i][j] = 0;
   }

   for (j=0 ; j<gene_num ; j++){
     if (M_q[j]==1){
       cluster_go_count      [E_q[j]]++;
       gocat_cluster_count[i][E_q[j]]++;     
       cat_gene_count++ ;
     }
   }

   free(M_q);

   for (j=0 ; j<cluster_num ; j++){

    //printf("ov=%d, s1=%d, s2=%d, N=%d\n", cluster_go_count[j], cluster_gene_count[j], cat_gene_count, gene_num);

    pvalue  [j][i] = log10(cumhyper  (cluster_go_count[j], cluster_gene_count[j], cat_gene_count, gene_num));

    // printf("Done\n"); 
    pvalue_u[j][i] = log10(cumhyper_u(cluster_go_count[j], cluster_gene_count[j], cat_gene_count, gene_num));
   }

   //printf("P-values for cat %d calculated\n", i);

 }


 printf("Done\n") ;

 //
 //writing the p-value matrix to file
 //


 printf("Writing the p-value and count matrices ... ") ;

 f = fopen(pvaluematrixfile, "w");
 if (f == 0)
   die("Cannot open pvaluematrixfile\n");

 fc = fopen(countmatrixfile, "w");
 if (fc == 0)
   die("Cannot open countmatrixfile\n");

 fprintf(f,  "GO");
 fprintf(fc, "GO");

 // row of cluster indices
 for (i=0 ; i<cluster_num ; i++){
   if (quantized == 1) {
     fprintf(f,  "\t%d", i) ;
     fprintf(fc, "\t%d", i) ;
   } else {
     fprintf(f,  "\t[%3.2f %3.2f]", E_q_bins[i], E_q_bins[i+1]);
     fprintf(fc, "\t[%3.2f %3.2f]", E_q_bins[i], E_q_bins[i+1]);
   }
 }
 fprintf(f,  "\n");
 fprintf(fc, "\n");

 // go thru all retained categories
 for(j=0 ; j<nb_cat ; j++) {
   // print desc
   fprintf(f,  "%s %s", go_terms_list[ A_go_indices[j] ],  go_desc[A_go_indices[j]]) ;
   fprintf(fc, "%s %s", go_terms_list[ A_go_indices[j] ],  go_desc[A_go_indices[j]]) ;
   // print row of p-values

   int totalcatcount = 0;
   for (i=0 ; i<cluster_num ; i++){
     totalcatcount += gocat_cluster_count[j][i];
   }

   //
   for (i=0 ; i<cluster_num ; i++){
     fprintf(f,  "\t%4.3f/%4.3f", pvalue[i][j], pvalue_u[i][j]) ;
     //int expcatnum = (int)(0.5 + totalcatcount * ( cluster_gene_count[i] / (double)gene_num ) );
     double fexpcatnum = totalcatcount * ( cluster_gene_count[i] / (double)gene_num ) ;
     
     fprintf(fc, "\t%d(%5.4f)", gocat_cluster_count[j][i], fexpcatnum) ;
   }

   fprintf(f,  "\n") ;
   fprintf(fc, "\n") ;
 }

 fclose(f) ;
 fclose(fc) ;


 if (outprofiles != 0) {
   
   FILE* fpop = fopen(outprofiles, "w");
   if (fpop != 0) {
     
     // printf("Outputing profiles of signif pathways to ... %s \n", outprofiles);

     fprintf(fpop, "GENE");
     for(j=0 ; j<nb_cat ; j++) {
       fprintf(fpop, "\t%s %s", go_terms_list[ A_go_indices[j] ],  go_desc[A_go_indices[j]]) ;
     }

     fprintf(fpop, "\n");


     for (i=0; i<gene_num; i++) {       

       fprintf(fpop, "%s", rownames[i]);
       for(j=0 ; j<nb_cat ; j++) {	 
	 fprintf(fpop, "\t%d", get_entry_in_binarized_array( go_profile[i], A_go_indices[j]));
       }      
       fprintf(fpop, "\n");

     }


   } else {
     printf("Cannot open profiles file\n");
   }

 }

 printf("Done\n") ;


 if (logfile != 0) {
   fclose (fplog);
 }

 return 0;
}


//
//
//
double evalSeed(int ii, unsigned char** gene_kmer, int nborfs, double score, int mbins, int* E_q, int ebins, int shuffle, int jn, int jn_f, int jn_t)
{

 double    pass1, pass2;
 int*   M_q;


 M_q = c_matrix_column_binarized(gene_kmer, ii, nborfs);


 pass1 = my_max_rank_test(score, M_q, mbins, E_q, ebins, nborfs, shuffle) ;
 //printf("pass1=%d\n", pass1);
 if (pass1 > 1){
   free(M_q);
   return 1;
 }

 // at this stage, if jn == 0, seeds passes.
 if (jn == 0) {
   free(M_q);
   return pass1;
 }

 //printf("jn=%d\tjn_t=%d\n", jn, jn_t);
 if (jn_t != 0)
   pass2 = jacknife_max_rank_test(M_q, mbins, E_q, ebins, nborfs, shuffle, jn, jn_f, jn_t, 0, 0);
 else
   pass2 =1 ;
 //printf("pass2=%d\n", pass2);
 if (pass2 == 0) {
   free(M_q);
   return 1;
 }

 free(M_q);

 return pass1;
}



//
//  returns whether new motif B_q (B) correlates with the A motifs
//
double minCondInfoNormalized(int *A_go_indices, unsigned char** go_profile, int DA, int nA, int* B_q, int DB, int* E_q, int DE, int n, int shuffle, double minr, int v, int* midx)
{
  int    i, j;
  int*   AB_q;
  int    DAE, DAB;
  double mi_ab_e, mi_a_e, mi_a_b;
  double cmi;
  double minratio = 1e6;
  double ratio;
  int*   M_q;

 DAE    = DA * DE;
 DAB    = DA * DB;
 *midx  = -1;

 for (i=0; i<nA; i++) {  // nA is the number of already kept GO cats

   j           = A_go_indices[i];  // get index of existing category
   M_q         = c_matrix_column_binarized(go_profile, j, n); // get corresponding profile for category


   // calculate conditional information I(Cnew;E|Cexist)  = I(A;E|B)  =  I(A,B;E) - I(A;E)
   AB_q        = combineQuantizedVectors(M_q, B_q, n, DA, DB);
   mi_ab_e     = CalculateMIbasic       (AB_q,   E_q, n, DAB, DE);
   mi_a_e      = CalculateMIbasic       (M_q, E_q, n, DA, DE);
   cmi         = mi_ab_e - mi_a_e;

   // calculate MI I(Cnew;Cexist)
   mi_a_b      = CalculateMIbasic       (M_q, B_q, n, DA, DB);

   if (cmi < 1e-16)
     cmi = 0.0;
   else if (mi_a_b < 1e-16)
     mi_a_b = 1e-16;

   ratio       = cmi / mi_a_b;

   if (i == 0)
     minratio = ratio;
   else {

     if (ratio < minratio) {
       minratio = ratio;
     }
   }

   free(M_q);
   free(AB_q);

   //
   //  break early feature
   //
   if (minratio < minr) {
     *midx = i;  // index of category that killed the current one
     return minratio;
   }



 }

 return minratio;

}

void read_go_index_data(char *filename, int *gene_num, int **line_sizes, char**** go_terms, char*** gene_names, char*** go_cat_list, int go_terms_num, int cat_min_count, int cat_max_count, struct my_hsearch_data *hash_excluded)
{
 ENTRY e;
 ENTRY *ep;

 char* s ;
 int mym ;
 int myn ;
 char* buff ;
 FILE* fp ;
 int i,j ;
 int mynmax = 100000 ;

 int    *t_line_sizes ;
 char ***t_go_terms ;

 struct my_hsearch_data *hash_go_count ;

 int t_idx = 0 ;
 int *A_count ;
 int hashret = 0;

 printf("entering\n") ;

 hash_go_count = (struct my_hsearch_data*)calloc(1,  sizeof(struct my_hsearch_data));
 my_hcreate_r(100000, hash_go_count);

 buff = (char*)malloc(100000 * sizeof(char)) ;
 fp = fopen(filename, "r") ;
 if (!fp)
 {
   printf("read_go_index_data: please enter a valid filename (%s invalid)\n", filename) ;
   exit(0) ;
 }

 *gene_names = (char**) malloc(mynmax * sizeof(char*)) ;
 *go_terms   = (char***)malloc(mynmax * sizeof(char**)) ;
 t_go_terms   = (char***)malloc(mynmax * sizeof(char**)) ;
 *line_sizes = (int*)   malloc(mynmax * sizeof(int)) ;
 t_line_sizes = (int*)   malloc(mynmax * sizeof(int)) ;

 myn = 0 ;

 A_count = (int*)malloc(100000 * sizeof(int)) ;
 for (i=0 ; i<go_terms_num ; i++)
   A_count[i] = 0 ;

 while (!feof(fp))
 {

   fgets(buff, mynmax, fp);

   if (feof(fp))
     break;

   chomp(buff);

   //printf("buff = %s\n", buff);

   // estimate the number of columns in line
   mym = 0;

   s = mystrtok(buff, '\t'); free(s);
   while ((s = mystrtok(0, '\t')))
   {
     mym ++;
     free(s);
   }

   t_line_sizes[myn] = mym ;

   t_go_terms[myn] = (char**)malloc(mym * sizeof(char*)) ;

   // get gene name
   s = mystrtok(buff, '\t') ;
   (*gene_names)[myn] = strdup(s) ;

   free(s) ;

   // get GO terms
   i = 0;

   while ((s = mystrtok(0, '\t')))
   {
     t_go_terms[myn][i] = strdup(s) ;

     e.key = t_go_terms[myn][i] ;

     my_hsearch_r(e, FIND, &ep, hash_go_count);

     if (!ep){
       e.key = strdup(t_go_terms[myn][i]) ;
       e.data = (char*) t_idx ;
       //printf("passed: %s->%s", e.key, e.data) ;
       //getchar() ;


       hashret = my_hsearch_r(e, ENTER, &ep, hash_go_count);
       if (hashret == 0) {
	 printf("Could not create hash table ...\n");
	 exit(0);
       }
       t_idx++ ;
     }else{
       int index = (int) ep->data ;
       A_count[index]++ ;
     }

     free(s) ;
     i++ ;
   }

   myn++ ;

   if (myn == mynmax)
   {
     die("read_go_index_data: running out of memory, please recompile ..\n") ;
   }
 }

 *gene_num = myn ;
 fclose(fp) ;

 //printf("counting the genes in each go term:%d\n", cat_max_count) ;


 //
 // exclude pathways if they have too many or too little members
 //
 for (i=0 ; i<go_terms_num ; i++){
   e.key = strdup((*go_cat_list)[i]) ;
   my_hsearch_r(e, FIND, &ep, hash_go_count);

   if(!ep)
     continue ;
   int index = (int) ep->data ;
   int count = A_count[index] ;
   if ((count < cat_min_count) || (count>cat_max_count)) {
     e.key = strdup((*go_cat_list)[i]) ;
     e.data = (char*) 1 ;
     hashret = my_hsearch_r(e, ENTER, &ep, hash_excluded);
     if (hashret == 0) {
       printf("Could not create hash table ...\n");
       exit(0);
     }
   }
 }


 //
 // MAKE THE FINAL LIST OF PATHWAYS
 //
 for (i=0 ; i<*gene_num ; i++){
   (*line_sizes)[i] = 0 ;
   int size = t_line_sizes[i] ;
   (*go_terms)[i] = (char**)malloc(size * sizeof(char*)) ;
   for (j=0 ; j<t_line_sizes[i] ; j++){
     e.key = strdup(t_go_terms[i][j]) ;
     //printf("go term %s\n", t_go_terms[i][j]) ;
     my_hsearch_r(e, FIND, &ep, hash_excluded);

     if (!ep){
       //printf("cat %s not excluded\n", t_go_terms[i][j]) ;
       //getchar() ;
       (*go_terms)[i][(*line_sizes)[i]] = strdup(t_go_terms[i][j]) ;
       (*line_sizes)[i]++ ;
     }
     else{
       //printf("cat %s is excluded;\n",  e.key) ;
       //getchar() ;
       continue ;
     }
   }
 }

 my_hdestroy_r(hash_go_count);

}

void read_go_names_data(char *filename, int *go_terms_num, char*** go_terms, char*** go_desc, char* cat_status, struct my_hsearch_data *hash_excluded)
{
 char* s ;
 ENTRY e;
 ENTRY *ep;

 int myn ;
 char* buff ;
 FILE* fp ;
 int mynmax = 100000 ;
 char status ;
 int hashret = 0;
 printf("entering\n");

 buff = (char*)malloc(100000 * sizeof(char)) ;
 fp = fopen(filename, "r") ;
 if (!fp)
 {
   printf("read_go_names_data: please enter a valid filename (%s invalid)\n", filename) ;
   exit(0) ;
 }

 *go_desc =  (char**)malloc(mynmax * sizeof(char*)) ;
 *go_terms = (char**)malloc(mynmax * sizeof(char*)) ;

 myn = 0 ;
 while (!feof(fp))
 {

   fgets(buff, mynmax, fp);

   if (feof(fp)) {
     break;
   }

   chomp(buff);

   // get GO term
   s = mystrtok(buff, '\t') ;
   //printf("%s\n", s) ;
   (*go_terms)[myn] = strdup(s) ;
   //(*go_terms)[myn][10] = '\0' ;
   free(s) ;
   //printf("l = %d\n", strlen((*go_terms)[myn]));

   // get GO description
   s = mystrtok(0, '\t') ;
   //printf("%s\n", s) ;
   (*go_desc)[myn] = strdup(s) ;
   free(s) ;

   //get cat status
   s = mystrtok(0, '\t') ;
   if (s == 0) {
     s = strdup("P");
   }
   status = s[0] ;
   //printf("%c: %d -> %d\n", status, get_cat_id(status), cat_status[get_cat_id(status)]) ;

   free(s) ;

   // exclude if not the right type of pathway
   if (cat_status[get_cat_id(status)]==0){
     e.key = strdup((*go_terms)[myn]) ;
     e.data = (char*) 1 ;
     hashret = my_hsearch_r(e, ENTER, &ep, hash_excluded);
     if (hashret == 0) {
       printf("Could not create hash table ...\n");
       exit(0);
     }
     //printf("cat %s excluded.\n", e.key) ;
     //getchar() ;
   }

   //printf("%s\t%s\n", (*go_terms)[myn], (*go_desc)[myn]);

   myn++ ;

   if (myn == mynmax)
   {
     die("read_go_names_data: running out of memory, please recompile ..\n") ;
   }

 }

 *go_terms_num = myn ;
 //getchar() ;
 fclose(fp) ;
 //getchar() ;
}

int CmpFunc(const void* _a, const void* _b)
{
 const GO* a = (const GO*) _a;
 const GO* b = (const GO*) _b;

 if (a->score < b->score)
   return 1;
 else if(a->score == b->score)
   return  0;
 else
   return -1;
}

int get_cat_id(char c)
{
 if (c=='F' || c=='f')
   return 0 ;
 if (c=='P' || c=='p')
   return 1 ;
 if (c=='C' || c=='c')
   return 2 ;
 return -1 ;
}


double my_max_rank_test(double score, int* M_q, int mbins, int* E_q, int ebins, int nborfs, int shuffle)
{

  double val = 0 ;
  int      j;
  int*     E_q_shu;
  double    mymi;
  double*   a_mi_shu;
  double    c;

  a_mi_shu = (double*)malloc(shuffle * sizeof(double));

  for (j=0; j<shuffle; j++) {
    E_q_shu = shuffleInt(E_q, nborfs);
    mymi = CalculateMIbasic(M_q, E_q_shu, nborfs, mbins, ebins);
    a_mi_shu[j] = mymi;
    free(E_q_shu);
  }

  qsort((void*)a_mi_shu, shuffle, sizeof(double), CmpDblRegular);

  c  = a_mi_shu[ shuffle - 1 ];

  if ((c < 1e-10) && (score < 1e-10))  // shortcut
    val = shuffle;
  else if (score <= a_mi_shu[0]) {  // other shortcut
    val = shuffle;
  } else {
    j = shuffle - 1;
    while ((j >= 0) && (score <= a_mi_shu[j])) {   // go while the shuffled score is higher than the real one
      j --;
    }
    val = shuffle - j - 1;

  }

  free(a_mi_shu);

  //printf("%f\t%f\n", score, c);

  //if (score > c)
  //  return 1;
  //else
  //  return 0;
  val = val/shuffle ;
  return val;
}
