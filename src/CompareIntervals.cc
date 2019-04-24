#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <string>
#include <vector>
#include <dirent.h>

extern "C" {
#include "dataio.h"
#include "statistics.h"
#include "sequences.h"
}
#include "interval_tree.h"

using namespace std;

#define CHRY_INDX	CHR_NUM - 1	/* chrY index */
#define CHRX_INDX	CHR_NUM - 2	/* chrX index */
#define CHR_OFFSET	3			/* chromosome number offset */
#define NMAX		100000

class GenInt : public Interval
{
public:
	int i;			// start position
	int j;			// end position
	string c;		// chromosome
	string id;		// id
	float score;	// score
	int strand;		// strand
	string str;
	
	int GetLowPoint() const 
	{
		return i;
	}
	int GetHighPoint() const
	{
		return j;
	}
};

vector<GenInt> a_int;
int		peaks_ov_cnt; /* Stores the number of peaks from file1(and random files) that have overlaps */

int searchChromsArray(const char* chr, char** chromsArray, int numchroms);
void readIntervals(char *intervals, int iswig, int hasheader, int hasid,
				   int ext1, int ext2, int show_peakdesc,
				   int showscores, int showstrand, int *numint);

void compareIntervals(char **chrnames,
					  int numchroms, int usepeak1center, IntervalTree* chr_tree,
					  int showprofile, int showovpeaks, int show_peakdesc,
					  int show_peakdesc2, int showscores, int showstrand,
					  int show_ov_int, int showunion, int included, int ovtype, 
					  int hasid1, int hasid2, char *outfile, FILE *outputFile,
					  int print_random_peaks);

int main(int argc, char** argv) {
	
	// files to compare
	FILE* f2;	
	
	/* directory for the random intervals files */
	DIR             *dip = NULL;
	struct dirent   *dit = NULL;
	
	// reading file
	char* buff			= NULL;
	char** p			= NULL;
	int    m			= 0;
	
	char* intervals1	= NULL;
	char* intervals2	= NULL;
	char* rand_dir		= NULL;
	char* prefix		= NULL;
	char* outfile		= NULL;
	FILE* outputFile	= NULL;
	
	float	rand_ovs_array[NMAX];
	
	int   verbose = 0;
	
	// store intervals
	IntervalTree* chr_tree;
	vector<GenInt> a_int2;
	GenInt tmpInt;
	
	int     numint1 = 0;
	int     numint2 = 0;
	
	int     ext1   = 0;
	int     ext2   = 0;
	int     hasid1 = 0;
	int     hasid2 = 0;
	int     rand2  = 0;
	
	int     showscores	= 0;
	int     showprofile = 0;
	int     show_ov_int = 0;
	char*   featname	= NULL;
	
	int     iswig		= 0;
	int     showunion	= 0;
	int     hasheader1	= 0;
	int     hasheader2	= 0;
	
	int		showstrand			= 0; // print strand
	int     usepeak1center		= 0; // convert a peak-interval into its center ?
	int     usepeak2center		= 0; // convert a peak-interval into its center ?
	int     usemaxpeak2height	= 0;
	int     showovpeaks		= 0;
	int     ovtype			= 0;	// 0 = AND, 1 = AND NOT
	char*   ovtype_text		= NULL;
	char*   outformat		= NULL;
	int     show_peakdesc	= 0;
	int     show_peakdesc2	= 0;	// show peak description of second file
	int     included		= 0;
	
	int		print_random_peaks	= 1;	/* we set this variable to 0 before running compare for random intervals: don't print the overlaps with the random peaks (we only want the count) */
	int		ov_peaks_file1		= 0;	/* stores the overlaps count only from file1 */
	int		all_peaks_file1		= 0;	/* stores the number of peaks only from file1*/
	
	char**  chrnames        = NULL;
	int		numchroms		= 24;		//number of chromosomes
	char	path[PATH_MAX];
	int 	i				= 0;
	int		cnt				= 0;
	vector<GenInt*> a_int_tmp;
	GenInt* bar;
	
	if (argc < 2) {
		die((char*)"Usage: ComparePeaks -peakfile1 FILE -peakfile2 FILE  [ -ext1 INT -hasid1 INT -ext2 INT -hasid2 INT -show_ov_int INT -showprofile INT -ovtype and|andnot -output peaklist|profile ]\n");
	}
	
	/* Store input options */
	if (exist_parameter(argc, argv, (char*)"-intervals1"))
		intervals1 = get_parameter(argc, argv, (char*)"-intervals1");
	
	if (exist_parameter(argc, argv, (char*)"-peakfile1"))
		intervals1 = get_parameter(argc, argv, (char*)"-peakfile1");
	
	if (exist_parameter(argc, argv, (char*)"-ext1"))
		ext1 = atoi(get_parameter(argc, argv, (char*)"-ext1"));
	
	if (exist_parameter(argc, argv, (char*)"-hasid1"))
		hasid1 = atoi(get_parameter(argc, argv, (char*)"-hasid1"));
	
	if (exist_parameter(argc, argv, (char*)"-hasheader1"))
		hasheader1 = atoi(get_parameter(argc, argv, (char*)"-hasheader1"));
	
	if (exist_parameter(argc, argv, (char*)"-ext2"))
		ext2 = atoi(get_parameter(argc, argv, (char*)"-ext2"));
	
	if (exist_parameter(argc, argv, (char*)"-hasid2"))
		hasid2 = atoi(get_parameter(argc, argv, (char*)"-hasid2"));
	
	if (exist_parameter(argc, argv, (char*)"-hasheader2"))
		hasheader2 = atoi(get_parameter(argc, argv, (char*)"-hasheader2"));
	
	if (exist_parameter(argc, argv, (char*)"-rand2"))
		rand2 = atoi(get_parameter(argc, argv, (char*)"-rand2"));
	
	if (exist_parameter(argc, argv, (char*)"-intervals2"))
		intervals2 = get_parameter(argc, argv, (char*)"-intervals2");
	
	if (exist_parameter(argc, argv, (char*)"-peakfile2"))
		intervals2 = get_parameter(argc, argv, (char*)"-peakfile2");
	
	if (exist_parameter(argc, argv, (char*)"-verbose"))
		verbose = atoi(get_parameter(argc, argv, (char*)"-verbose"));
	
	if (exist_parameter(argc, argv, (char*)"-showscores"))
		showscores = atoi(get_parameter(argc, argv, (char*)"-showscores"));
	
	if (exist_parameter(argc, argv, (char*)"-showprofile"))
		showprofile = atoi(get_parameter(argc, argv, (char*)"-showprofile"));
	
	if (exist_parameter(argc, argv, (char*)"-showovpeaks"))
		showovpeaks = atoi(get_parameter(argc, argv, (char*)"-showovpeaks"));
	
	if (exist_parameter(argc, argv, (char*)"-show_ov_int"))
		show_ov_int = atoi(get_parameter(argc, argv, (char*)"-show_ov_int"));
	
	if (exist_parameter(argc, argv, (char*)"-2includes1"))
		included = atoi(get_parameter(argc, argv, (char*)"-2includes1"));
	
	if (exist_parameter(argc, argv, (char*)"-featname"))
		featname = get_parameter(argc, argv, (char*)"-featname");
	
	if (exist_parameter(argc, argv, (char*)"-output")) {
		outformat = get_parameter(argc, argv, (char*)"-output");
		if (strcmp(outformat, "peaklist") == 0) {
			showovpeaks = 1;
		} else if (strcmp(outformat, "profile") == 0) {
			showprofile = 1;
		}	
	}
	
	if (exist_parameter(argc, argv, (char*)"-iswig"))
		iswig = atoi(get_parameter(argc, argv, (char*)"-iswig"));
	
	if (exist_parameter(argc, argv, (char*)"-showunion"))
		showunion = atoi(get_parameter(argc, argv, (char*)"-showunion"));
	
	if (exist_parameter(argc, argv, (char*)"-showstrand"))
		showstrand = atoi(get_parameter(argc, argv, (char*)"-showstrand"));
	
	if (exist_parameter(argc, argv, (char*)"-usepeak1center"))
		usepeak1center = atoi(get_parameter(argc, argv, (char*)"-usepeak1center"));
	
	if (exist_parameter(argc, argv, (char*)"-usepeak2center"))
		usepeak2center = atoi(get_parameter(argc, argv, (char*)"-usepeak2center"));
	
	if (exist_parameter(argc, argv, (char*)"-usemaxpeak2height"))
		usemaxpeak2height = atoi(get_parameter(argc, argv, (char*)"-usemaxpeak2height"));
	
	if (exist_parameter(argc, argv, (char*)"-showpeakdesc"))
		show_peakdesc = atoi(get_parameter(argc, argv, (char*)"-showpeakdesc"));
	
	if (exist_parameter(argc, argv, (char*)"-showpeakdesc2"))
		show_peakdesc2 = atoi(get_parameter(argc, argv, (char*)"-showpeakdesc2"));
	
	
	if (exist_parameter(argc, argv, (char*)"-ovtype")) {
		ovtype_text = uc(get_parameter(argc, argv, (char*)"-ovtype"));
		if (strcmp(ovtype_text, "AND") == 0) {
			ovtype = 0;
		} else if (strcmp(ovtype_text, "ANDNOT") == 0) {
			ovtype = 1;
		} else {
			die((char*)"ovtype unknown.\n");
		}
	}
	
	if (exist_parameter(argc, argv, (char*)"-outfile")) {
		outfile = get_parameter(argc, argv, (char*)"-outfile");
	}
	
	if (exist_parameter(argc, argv, (char*)"-rand_dir"))
		rand_dir = get_parameter(argc, argv, (char*)"-rand_dir");
	
	if (exist_parameter(argc, argv, (char*)"-prefix"))
		prefix = get_parameter(argc, argv, (char*)"-prefix");	
	
	if(outfile != 0) {
		outputFile = fopen(outfile, "w");
	}
	
	/*if (outputFile==NULL)
	{
		fprintf(stdout, "Could not open outfile: %s. Please try again.\n", outfile);
		exit(1);
	}*/
	
	readChrData2(intervals1, intervals2, hasid1, hasid2, &chrnames, &numchroms);
	
	if (verbose == 1) {
		(void)fprintf(stdout, "Found %d chroms\n", numchroms);
	}
	
	/* Initialize the Interval trees array */
	chr_tree = new IntervalTree[numchroms];
	
	default_set_seed(time(NULL));
	
	/* Read the first file of intervals. Store intervals in vector a_int. */
	/* This needed to be a function so that we can call it many times, one for each random intervals file. */
	readIntervals(intervals1, iswig, hasheader1, hasid1, ext1, ext2, show_peakdesc, 
				  showscores, showstrand, &numint1);
	
	/* 
	 *
	 * Read the second file of intervals and build one interval tree per chromosome
	 *
	 */
	
	numint2 = 0;
	
	/* Open second file */
	f2 = fopen(intervals2, "r");
	
	if (!f2) {
		(void)fprintf(stdout, "Cannot open peakfile2 (%s)\n", intervals2);
		exit(1);
	}
	
	if ((iswig == 1) || (hasheader2 == 1))
		fgets(buff, NMAX, f2);
	
	buff  = (char*)malloc(NMAX * sizeof(char));  
	
	/* read intervals line by line */
	while (!feof(f2)) {
		fgets(buff, NMAX, f2);
		if (feof(f2))
			break; 
		chomp(buff);
		char* line = strdup(buff);
		
		split_line_delim(buff, (char*)"\t", &p, &m);
		
		// get a new entry
		a_int2.push_back(tmpInt);
		
		/* store interval attributes in vector a_int2 */
		int idx = 0;
		if (hasid2 == 1) {
			idx = 1;
			a_int2[numint2].id = string(p[0]);
		}
		
		a_int2[numint2].c = string(p[idx]);
		
		int p1 = atoi(p[idx+1]) - ext2;
		int p2 = atoi(p[idx+2]) + ext2;
		
		// transform interval into single point
		if (usepeak2center == 1) {
			p1 = (p1 + p2)/2;
			p2 = p1;
		} else if (usemaxpeak2height == 1) {
			if (idx + 5 < m) {
				int pkc = atoi(p[idx+5]);
				if ((p1 <= pkc) && (pkc <= p2)) {
					p1 = pkc;
					p2 = pkc;
				} else {
					die((char*)"max peak info not between peak boundaries, abort. Please verify -peak2file.\n");
				}
			} else {
				die((char*)"max peak info not present in -peak2file.\n");
			}
			
		}
		
		/* if rand2 is 1, we create a random interval */
		if (rand2 == 1) {
			int li = p2 - p1 + 1;
			p1 = rand_chr_interval(p[idx], li);
			p2 = p1 + li;
		}
		
		/* store interval attributes in vector a_int2 */
		a_int2[numint2].i = p1;
		a_int2[numint2].j = p2;
		
		if (show_peakdesc == 1) {
			a_int2[numint2].str = string(line);
		} else {
			free(line);
		}
		
		if (showscores == 1) {
			a_int2[numint2].score = atof(p[idx+4]);
		}
		
		/* Find t_indx, which corresponds to the chromosome of the current interval */
		int t_indx  = 0;
		
		t_indx = searchChromsArray(a_int2[numint2].c.c_str(), chrnames, numchroms);
		
		if (t_indx == -1) {
			continue;
		}
		
		/* Insert interval in the corresponding interval tree */
		bar = new GenInt(a_int2[numint2]);
		a_int_tmp.push_back(bar);
		chr_tree[t_indx].Insert(bar);
		
		free(p);
		numint2++;
	}
	
	/* Close second file */
	(void)fclose(f2);
	free(buff);
	
	if (numint2 == 0) {
		(void)fprintf(stdout, "No peaks found in %s\n", intervals2);
	}
	if (verbose == 1) {
		(void)fprintf(stdout, "Both files were read.\n");
	}
	
	if ((showprofile >= 1) && (featname != 0)) {
		if(outfile != 0) {
			(void)fprintf(outputFile, "FEATURE\t%s\nFEATTYPES\t%d\n", featname, (showprofile==1?0:1));
		}
		else {
			printf("FEATURE\t%s\nFEATTYPES\t%d\n", featname, (showprofile==1?0:1));
		}
	}
	
	if (showprofile == 1) 
		if(outfile != 0) {
			(void)fprintf(outputFile, "GENE\tBOUND\n");
		}
		else {
			printf("GENE\tBOUND\n");
		}
	
	peaks_ov_cnt = 0; /* it will get a value from the function compareIntervals */
	
	/* This function searches the interval trees for the intervals of file1 (stored in vector). */
	/* This needed to be a function so that we can call it many times, one for each random intervals file.*/
	compareIntervals(chrnames, numchroms, usepeak1center,
					 chr_tree, showprofile, showovpeaks, show_peakdesc,
					 show_peakdesc2, showscores, showstrand, show_ov_int, showunion,
					 included, ovtype, hasid1, hasid2, outfile, outputFile, print_random_peaks);
	
	all_peaks_file1	= numint1;

	if (ovtype == 0) {
		ov_peaks_file1	= peaks_ov_cnt;
	}
	else {
		ov_peaks_file1	= all_peaks_file1 - peaks_ov_cnt;
	}
		
	/* If rand_dir is set, open directory and read random intervals files */
	if (rand_dir != NULL && prefix != NULL) {
		
		if ((dip = opendir(rand_dir)) == NULL)
		{
			(void)fprintf(stdout, "Cannot open rand_dir (%s)\n", rand_dir);
			exit(1);
		}
		
		int i = 0; // this variable is used to store the peaks_ov_count in array
		
		/* For each random intervals file, read it (store in a_int), and look for overlaps with the interval tree. */
		while ((dit = readdir(dip)) != NULL)
		{
			if (strncmp(dit->d_name, prefix, strlen(prefix)) == 0) {
				
				/* Reset vector and numint counter. */
				a_int.clear();
				numint1 = 0;
				
				/* generate the full path. */
				(void)memset(path, 0, PATH_MAX);
				(void)memcpy(path, rand_dir, strlen(rand_dir));
				(void)strcat(path, "/");
				(void)strcat(path, dit->d_name);
				
				/* Read the file and store in a_int. */
				readIntervals(path, iswig, hasheader1, hasid1, ext1, ext2, show_peakdesc, 
							  showscores, showstrand, &numint1);
				
				peaks_ov_cnt		= 0;
				print_random_peaks	= 0;
				
				/* Look for overlaps of a_int with the intervals tree. */
				compareIntervals(chrnames, numchroms, usepeak1center,
								 chr_tree, showprofile, showovpeaks, show_peakdesc,
								 show_peakdesc2, showscores, showstrand, show_ov_int, showunion,
								 included, ovtype, hasid1, hasid2, outfile, outputFile, print_random_peaks);
							
				/* Store the overlaps count in array */
				
				if (ovtype == 0) {
					rand_ovs_array[i] = peaks_ov_cnt;
					if (rand_ovs_array[i] >= ov_peaks_file1)
						cnt ++;
				}
				else {
					rand_ovs_array[i] = numint1 - peaks_ov_cnt;
				}
												
				/* Increase array index */
				i++;
			}
		}
		
		/* Estimate  */
		float avg = average(rand_ovs_array, i);
		float std = stddev(rand_ovs_array, i); 
				
		float zscore = (ov_peaks_file1 - avg)/std;
		//float tscore = 10*zscore + 50;
		
		float pvalue = (float)cnt/i;

		(void)fprintf(stdout, "Observed: %d\t", ov_peaks_file1);
		(void)fprintf(stdout, "Average: %f\t", avg);
		(void)fprintf(stdout, "St. Deviation: %f\t", std);
		(void)fprintf(stdout, "Z-score: %f\t", zscore);	
		//(void)fprintf(stdout, "T-score: %f\n", tscore);	
		(void)fprintf(stdout, "P-value: %f\n", pvalue);
		
		/* Close directory */
		(void)closedir(dip);
	}
	
	/* delete tree */
	delete[] chr_tree;
	
	/* cleanup */
	for (i=0; i<numchroms; i++)
		free(chrnames[i]);

	for (i=0; i<a_int_tmp.size(); i++)
		delete a_int_tmp[i];

	free(chrnames);
	free(ovtype_text);
	
	return 0;
}

void compareIntervals(char **chrnames,
					  int numchroms, int usepeak1center, IntervalTree* chr_tree,
					  int showprofile, int showovpeaks, int show_peakdesc,
					  int show_peakdesc2, int showscores, int showstrand,
					  int show_ov_int, int showunion, int included, int ovtype, 
					  int hasid1, int hasid2, char *outfile, FILE *outputFile,
					  int print_random_peaks)
{
	
	int     x1, x2;
	int     minx;
	int     maxx;
	
	/* For each interval stored in vector a_int (intervals of first file) */
	for (size_t i=0; i < a_int.size(); i++) {
		
		minx = a_int[i].i;
		maxx = a_int[i].j;
		
		/* Find t_indx, which corresponds to the chromosome of the current interval */
		int t_indx	= 0;
		
		t_indx = searchChromsArray(a_int[i].c.c_str(), chrnames, numchroms);
		
		if (t_indx == -1) { // || a_int[i].c.find('_') != -1) {
			continue;
		}
		
		// use for overlap
		x1 = a_int[i].i;
		x2 = a_int[i].j;
		
		// transform interval into single point
		if (usepeak1center == 1) {
			x1 = (x1+x2)/2;
			x2 = x1;
		}
		
		/* Enumerate function returns all intervals stored in the tree that overlap with x1, x2 */
		TemplateStack<GenInt*> *res = (TemplateStack<GenInt*>*)chr_tree[t_indx].Enumerate(x1, x2);
		
		/* Number of overlapping intervals */
		int mycnt		= res->Size();
		
		if (mycnt != 0) {
			peaks_ov_cnt++;
		}
		
		// show profile type of output
		if ((showprofile == 1) && (mycnt >= 1))
			mycnt = 1;
		
		if (print_random_peaks !=0 ) {
			
			
			// show data for all peaks
			if (showovpeaks == 0) {
				if (show_peakdesc == 0) { 
					if (hasid1 == 1) {
						if(outfile != 0) {
							(void)fprintf(outputFile, "%s\t", a_int[i].id.c_str());
						}
						else {
							printf("%s\t", a_int[i].id.c_str());
						}
					}
					if (showprofile == 0) {
						if(outfile != 0) {
							(void)fprintf(outputFile, "%s\t%d\t%d\t", a_int[i].c.c_str(), a_int[i].i, a_int[i].j);
						}
						else {
							printf("%s\t%d\t%d\t", a_int[i].c.c_str(), a_int[i].i, a_int[i].j);
						}
					}
					if (showscores == 1) {
						if(outfile != 0) {
							(void)fprintf(outputFile, "%f\t", a_int[i].score);
						}
						else {
							printf("%f\t", a_int[i].score);
						}
					}
					if (showstrand == 1) {
						if(outfile != 0) {
							(void)fprintf(outputFile, "%d\t", a_int[i].strand);
						}
						else {
							printf("%d\t", a_int[i].strand);
						}
					}
				}
				else { // show peak info
					if(outfile != 0) {
						(void)fprintf(outputFile, "%s\t", a_int[i].str.c_str());
					}
					else {
						printf("%s\t", a_int[i].str.c_str());
					}
				}
				
				//print number of overlaps
				if(outfile != 0) {
					(void)fprintf(outputFile, "%d", mycnt);
				}
				else {
					printf("%d", mycnt);
				}
				
				// show the overlapping peaks
				if (show_ov_int == 1) {
					for (int it = 0; it < res->Size(); it++) {
						if ((*res)[it]->i < minx)
							minx = (*res)[it]->i;
						if ((*res)[it]->j > maxx)
							maxx = (*res)[it]->j;
						
						if(outfile != 0) {
							(void)fprintf(outputFile, "\t%s-%d-%d", (*res)[it]->c.c_str(), (*res)[it]->i, (*res)[it]->j);
						}
						else {
							printf("\t%s-%d-%d", (*res)[it]->c.c_str(), (*res)[it]->i, (*res)[it]->j);
						}					
						
						if (hasid2 == 1) {
							if(outfile != 0) {
								(void)fprintf(outputFile, ":%s", (*res)[it]->id.c_str());
							}
							else {
								printf(":%s", (*res)[it]->id.c_str());
							}
						}
					}
				}
				
				if (showunion == 1) {
					if(outfile != 0) {
						(void)fprintf(outputFile, "\t%s-%d-%d", a_int[i].c.c_str(), minx, maxx);
					}
					else {
						printf("\t%s-%d-%d", a_int[i].c.c_str(), minx, maxx);
					}
				}
				
				if(outfile != 0) {
					(void)fprintf(outputFile, "\n");
				}
				else {
					printf("\n");
				}
				
			} 
			// show only peaks that overlap (showovpeaks == 1)
			else {
				// refine overlap data if inclusion required
				if (included == 1) {
					int newcnt = 0; 
					for (int it = 0; it < res->Size(); it++) {
						// interval includes peak
						if (((*res)[it]->i <= a_int[i].i) && ( a_int[i].j <= (*res)[it]->j )) {
							newcnt ++;
						}
					}
					mycnt = newcnt;
				}
				
				if ((ovtype == 0) && (mycnt > 0)) {  // AND
					
					if (show_peakdesc == 0) {  // reduce peak rep
						if (hasid1 == 1) {
							if(outfile != 0) {
								(void)fprintf(outputFile, "%s\t", a_int[i].id.c_str());
							}
							else {
								printf("%s\t", a_int[i].id.c_str());
							}
						}
						
						if(outfile != 0) {
							(void)fprintf(outputFile, "%s\t%d\t%d", a_int[i].c.c_str(), a_int[i].i, a_int[i].j);
						}
						else {
							printf("%s\t%d\t%d", a_int[i].c.c_str(), a_int[i].i, a_int[i].j);				    
						}
						
						if (showscores == 1) {
							if(outfile != 0) {
								(void)fprintf(outputFile, "\t%f\t", a_int[i].score);
							}
							else {
								printf("\t%f\t", a_int[i].score);
							}
						}
						
					} else { // show peak info
						if(outfile != 0) {
							(void)fprintf(outputFile, "%s", a_int[i].str.c_str());
						}
						else {
							printf("%s", a_int[i].str.c_str());
						}
					}
					
					//print number of overlaps
					
					/*if(outfile != 0) {
					 (void)fprintf(outputFile, "\t%d", mycnt);
					 }
					 else {
					 printf("\t%d", mycnt);
					 }*/
					
					
					// show the overlapping peaks
					if (show_ov_int == 1) {
						for (int it = 0; it < res->Size(); it++) {
							if (show_peakdesc2 == 0) {  // reduce peak rep
								if ((*res)[it]->i < minx)
									minx = (*res)[it]->i;
								if ((*res)[it]->j > maxx)
									maxx = (*res)[it]->j;
								
								if(outfile != 0) {
									(void)fprintf(outputFile, "\t%s-%d-%d", (*res)[it]->c.c_str(), (*res)[it]->i, (*res)[it]->j);
								}
								else {
									printf("\t%s-%d-%d", (*res)[it]->c.c_str(), (*res)[it]->i, (*res)[it]->j);
								}
								
							} else { // show peak info
								if(outfile != 0) {
									(void)fprintf(outputFile, "\t%s", (*res)[it]->str.c_str());
								}
								else {
									printf("\t%s", (*res)[it]->str.c_str());
								}
							}
							
							if (hasid2 == 1) {
								if(outfile != 0) {
									(void)fprintf(outputFile, ":%s", (*res)[it]->id.c_str());
								}
								else {
									printf(":%s", (*res)[it]->id.c_str());
								}
							}
							
							if (showscores == 1) {
								if(outfile != 0) {
									(void)fprintf(outputFile, ":%f", (*res)[it]->score);
								}
								else {
									printf(":%f", (*res)[it]->score);
								}
							}
						}
					}
					
					if (showunion == 1) {
						if(outfile != 0) {
							(void)fprintf(outputFile, "\t%s-%d-%d", a_int[i].c.c_str(), minx, maxx);
						}
						else {
							printf("\t%s-%d-%d", a_int[i].c.c_str(), minx, maxx);
						}
					}
					
					if(outfile != 0) {
						(void)fprintf(outputFile, "\n");
					}
					else {
						printf("\n");
					}
					
				} else if ((ovtype == 1) && (mycnt == 0)) {  // AND NOT
					
					if (show_peakdesc == 0) {  // reduce peak rep
						if(outfile != 0) {
							(void)fprintf(outputFile, "%s\t%d\t%d", a_int[i].c.c_str(), a_int[i].i, a_int[i].j);
						}
						else {
							printf("%s\t%d\t%d", a_int[i].c.c_str(), a_int[i].i, a_int[i].j);				    
						}
					} else {
						if(outfile != 0) {
							(void)fprintf(outputFile, "%s", a_int[i].str.c_str());
						}
						else {
							printf("%s", a_int[i].str.c_str());				      
						}
					}	
					
					if(outfile != 0) {
						(void)fprintf(outputFile, "\n");
					}
					else {
						printf("\n");
					}			
				}
			}
		}

		/* cleanup */
		delete res;
	}
}

void readIntervals(char *intervals, int iswig, int hasheader, int hasid,
				   int ext1, int ext2, int show_peakdesc, int showscores,
				   int showstrand, int *numint)
{
	FILE *f1	= NULL;
	char *buff	= NULL;
	char** p	= NULL;
	int    m	= 0;
	GenInt tmpInt;
	
	buff  = (char*)malloc(NMAX * sizeof(char));  
	
	// read first set of intervals
	*numint = 0;
	
	/* Open first file	*/
	f1 = fopen(intervals, "r");
	if (!f1) {
		(void)fprintf(stdout, "Cannot open peakfile1 (%s)\n", intervals);
		exit(1);
	}
	
	if ((iswig == 1) || (hasheader == 1))
		fgets(buff, NMAX, f1);
	
	/* read intervals line by line */
	while (fgets(buff, NMAX, f1) != 0) {
		
		chomp(buff);
		char* line = strdup(buff);
		
		split_line_delim(buff, (char*)"\t", &p, &m);
		
		// get a new entry
		a_int.push_back(tmpInt);
		
		/* store interval attributes in vector a_int2 */
		int idx = 0;
		if (hasid == 1) {
			idx = 1;
			a_int[*numint].id = string(p[0]);
		} 
		a_int[*numint].c = string(p[idx]);
		a_int[*numint].i = atoi(p[idx+1]) - ext1;
		a_int[*numint].j = atoi(p[idx+2]) + ext1;
		
		if (show_peakdesc == 1) {
			a_int[*numint].str = string(line);
		} else {
			free(line); /// no need anymore
		}
		
		if (showscores == 1) {
			a_int[*numint].score = atof(p[idx+4]);
		}
		if (showstrand == 1) {
			a_int[*numint].strand = atoi(p[idx+3]);
		}
		free(p);
		(*numint)++;
	}
	
	if (*numint == 0) {
		(void)fprintf(stdout, "No peaks found in %s\n", intervals);
	}
	
	/* Close first file */
	(void)fclose(f1);
	free(buff);
}

int searchChromsArray(const char* chr, char** chromsArray, int numchroms)
{	
	for (int i=0; i<numchroms; i++) {
		if (strcmp(chr, chromsArray[i]) == 0) {
			return i;
		}
	}
	return -1;
}
