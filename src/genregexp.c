#include <stdio.h>

#ifdef BNS
#include <pcre/pcre.h>
#else
#include <pcre.h>
#endif

#include <string.h>
#include <ctype.h>

#define max(x,y) ((x)>(y)?(x):(y))
#define min(x,y) ((x)<(y)?(x):(y))

//functions
int		findSites(char* r, char* name, char* seq, int* positions, int* np, int* orientations, int* no, int* lengths_m, int singlestrand, int lenoffset); 
char	*substr(char *string, int start, int length); 
void	chomp(char* s); 
char*	complement(char* s); 
char*	get_parameter(int argc, char** argv, char* param);
char*	nextSequence(FILE* fp, char** name, int* size); 

//variables
int		verbose						= 1;
char*	nextSequence_currentLine; 
int		nextSequence_started		= 0;
int		nextSequence_ended			= 0;
int		recountpalindromes;
int		reldist						= 0;

//main function
int main(int argc, char** argv) {
	
	int	i	= 0;
	int	l	= 0;
	FILE* fp;
	
	char* seq	= NULL;
	char* name	= NULL;
	int   size	= 0;
	
	int*  positions		= NULL;
	int*  orientations	= NULL;
	int*  lengths_m		= NULL;
	
	int   np	= 0;
	int   no	= 0;
	int   singlestrand	= 0;
	char* fastafile		= NULL;
	char* re			= NULL;
	int   simple		= 0;
	int   onlymatches	= 1;
	int   lenoffset		= 0;
	char* ss			= NULL;
	int   j				= 0;
	int   sp			= 0;
	int	  sl			= 0;;
	
	int   flank			= 0;
	int   flank_u		= 0;
	int   flank_d		= 0;
	int   po			= 0;
	double prct			= 0.0;
	int   format		= 0;
	char* formattext	= 0;
	
	positions    = (int*)malloc(5000 * sizeof(char));
	orientations = (int*)malloc(5000 * sizeof(char));
	lengths_m    = (int*)malloc(5000 * sizeof(char));
	
	nextSequence_currentLine = (char*)malloc(50000 * sizeof(char));
	
	//check for lack of arguments
	if (argc <= 1) {
		printf("Args: -re STR -fastafile FILE [ -rna INT ]\n");
		exit(1);
	}
	
	//
	// get parameters
	//
	fastafile  = get_parameter(argc, argv, "-fastafile");
	
	if (strlen(get_parameter(argc, argv, "-simple")) > 0) { 
		simple   = atoi(get_parameter(argc, argv, "-simple"));
		format   = 1;
	}
	
	if (strlen(get_parameter(argc, argv, "-format")) > 0) {
		formattext = get_parameter(argc, argv, "-format");    
		if (strcmp(formattext, "simple") == 0) {
			simple   = 1;
			format   = 1;
		} else if (strcmp(formattext, "count") == 0) {
			format   = 2;
		} else if (strcmp(formattext, "density") == 0) {
			format   = 3;
		}
	}
		
	if (strlen(get_parameter(argc, argv, "-lenoffset")) > 0) 
		lenoffset = atoi(get_parameter(argc, argv, "-lenoffset"));
	
	// takes priority
	if (strlen(get_parameter(argc, argv, "-flank")) > 0) {
		flank   = atoi(get_parameter(argc, argv, "-flank"));
		flank_u = flank;
		flank_d = flank;
	} else {
		
		if (strlen(get_parameter(argc, argv, "-flank_u")) > 0) 
			flank_u = atoi(get_parameter(argc, argv, "-flank_u"));
		
		if (strlen(get_parameter(argc, argv, "-flank_d")) > 0) 
			flank_d = atoi(get_parameter(argc, argv, "-flank_d"));
	}
	
	
	if (strlen(get_parameter(argc, argv, "-onlymatches")) > 0) 
		onlymatches   = atoi(get_parameter(argc, argv, "-onlymatches"));
	
	re = get_parameter(argc, argv, "-re");
	
	
	recountpalindromes = 0;
	if (strlen(get_parameter(argc, argv, "-recountpalindromes")) > 0) {
		recountpalindromes = atoi(get_parameter(argc, argv, "-recountpalindromes"));
	} else {
		recountpalindromes = 0;
	}
	
	//  relative distance
	if (strlen(get_parameter(argc, argv, "-reldist")) > 0) {
		reldist = atoi(get_parameter(argc, argv, "-reldist"));
	} else {
		reldist = 0;
	}
	
	if (strcmp(get_parameter(argc, argv, "-singlestrand"), "") != 0)
		singlestrand      = atoi(get_parameter(argc, argv, "-singlestrand"));
	else 
		singlestrand      = 0;
	
	if (strcmp(get_parameter(argc, argv, "-rna"), "") != 0)
		singlestrand      = atoi(get_parameter(argc, argv, "-rna"));
	else 
		singlestrand      = 0;
	
	//open fasta file
	fp = fopen(fastafile, "r");
	
	if (!fp) {
		printf("cannot open %s ..\n", fastafile);
		exit(0);
	}
	
	while ((seq = nextSequence(fp, &name, &size))) {   
		
		if (seq && (strlen(seq) > 0)) {
			
			np = 0;
			no = 0;
			
			l  = strlen(seq);
			//for (i=0; i<l; i++) {
			//seq[i] = toupper(seq[i]);
			//}
			
			findSites(re, name, seq, positions, &np, orientations, &no, lengths_m, singlestrand, lenoffset);
			
			if (format == 1) {
				printf("%s\t%d\n", name, (np>0?1:0));
			} 
			else if (format == 2) {
				printf("%s\t%d\n", name, np);
			} 
			else if (format == 3) {
				printf("%s\t%3.1f\n", name, 1000*np/(double)strlen(seq));
			} 
			else { // format == 0
				
				for (i=0; i<np; i++) {
					printf("%s\t", name);
					
					if (reldist == 1) {
						po = l - positions[i] - 1;
					} else {
						po = positions[i];
					}
					
					prct = (float)po/l*100;
					
					printf("%d\t%3.1f\t%d", po, prct, orientations[i]);
					
					fflush(stdout);
					
					if ((flank > 0) || (flank_u > 0) || (flank_d > 0)) {
						
						sp = max(0, positions[i] - flank_u);
						sl = positions[i] - sp;
						
						ss = substr(seq, sp, sl);
						for (j=0; j<flank_u; j++)
							ss[j] = tolower(ss[j]);
						printf("\t%s", ss);
						free(ss);
						
						fflush(stdout);
						
						ss = substr(seq, positions[i], lengths_m[i]);
						printf("%s", ss);
						free(ss);
						fflush(stdout);
						
						sp = positions[i] + lengths_m[i];
						
						//printf("\nsp=%d\n", sp);
						
						sl = min(l, sp+flank_d); 
						//printf("sl=%d, flank_d=%d\n", sl, flank_d);
						
						ss = substr(seq, sp, sl - sp);
						for (j=0; j<flank_d; j++)
							ss[j] = tolower(ss[j]);
						printf("%s", ss);
						free(ss);
						fflush(stdout);
					}
					
					printf("\n");
				}
			} 
		}
		free(seq);
		free(name);
		seq		= NULL;
		name	= NULL;
	}
	
	/* cleanup */
	fclose(fp);
	free(positions);
	free(orientations);
	free(lengths_m);
	free(nextSequence_currentLine);
	
	if(name != NULL)
		free(name);
	
	if(seq != NULL)
		free(seq);
	
	return 0;
}

/** REGEXP, returns 1 if a site was found, 0 otherwise **/
int findSites(char* r, char* name, char* seq, int* positions, int* np, int* orientations, int* no, int* lengths_m, int singlestrand, int lenoffset) 
{	
	pcre* re;
	const char* error;
	int   erroffset;
	int   rc;
	int   ovector[30];
	char* cr = 0;
	char*  substring_start;
	int  substring_length;
	int startoffset = 0;
	
	char* seq_pos;
	int pos;
	
	//
	//  store the position where a RE has been found, to avoid counting twice the palindromes
	//
	seq_pos = (char*)calloc(strlen(seq), sizeof(char));
	
	*np = 0;
	*no = 0;
	
	//printf("searching for %s in %s\n", r, seq);
	
	
	//
	// ONE STRAND
	//
	re = pcre_compile(r,
					  PCRE_CASELESS,
					  &error,
					  &erroffset,
					  NULL);
	
	while ( (rc = pcre_exec(re,
							NULL,
							seq,
							strlen(seq),
							startoffset,
							0,
							ovector,
							30) ) > 0) {
		
		
		substring_start    = seq         + ovector[0]; 
		substring_length   = ovector[1]  - ovector[0]; 
		
		// save the position
		pos = strlen(seq) - ovector[0] - 1;
		
		//if (reldist == 1) {
		//  positions[ *np ] = pos;
		//} else {
		positions[ *np ] = ovector[0];
		//}
		
		if (lengths_m != 0)
			lengths_m[ *np ] = substring_length;
		
		//printf("%s\n", substr(seq, pos, substring_length));
		//printf("%.*s\n", substring_length, substring_start);
		
		// save the orientation
		orientations[ *no ] = 1;
		
		if (lenoffset == 0)
			startoffset        = ovector[0] + substring_length;
		else
			startoffset        = ovector[0] + lenoffset;
		
		seq_pos[ pos ] = 1;
		
		(*np)++;
		(*no)++;
	}
	
	
	pcre_free(re);

	//return 0;
	
	if (singlestrand == 0) {
		
		// OTHER STRAND
		
		cr = complement(r);
		
		//printf("Complement of %s is %s\n", r, cr);
		
		startoffset = 0;
		
		re = pcre_compile(cr,
						  PCRE_CASELESS,
						  &error,
						  &erroffset,
						  NULL);
		
		while ( (rc = pcre_exec(re,
								NULL,
								seq,
								strlen(seq),
								startoffset,
								0,
								ovector,
								30)) > 0) {
			substring_start = seq + ovector[0]; 
			substring_length = ovector[1] - ovector[0]; 
			//printf("%.*s\n", substring_length, substring_start); 
			//  if not a palindromic version
			pos = strlen(seq) - ovector[0] - 1;
			if ((seq_pos[ pos ] != 1) || ((seq_pos[ pos ] == 1) && (recountpalindromes == 1)) ) {
				// save the position
				//if (reldist == 1) {
				//  positions[ *np ] = pos;
				// } else {
				positions[ *np ] = ovector[0];
				// }
				
				if (lengths_m != 0)
					lengths_m[ *np ] = substring_length;
				
				//positions[ *np ] = pos;
				(*np)++;
				
			} 
			
			// save the orientation in all cases
			orientations[ *no ] = -1;
			(*no)++;
			
			if (lenoffset == 0)
				startoffset        = ovector[0] + substring_length;
			else
				startoffset        = ovector[0] + lenoffset;
			
		}
		
	}
	
    
	if (cr != 0) 
		free(cr);
	
	free(seq_pos);
	pcre_free(re);
	return 0;
}

//
// take a char* and replace every nt by its complement
//
char* complement(char* s) 
{
	
	int l = strlen(s);
	char* c;
	char  d;
	int i;
	
	c = (char*)calloc(l+1, sizeof(char));
	
	for (i=l-1; i>=0; i--) {
		if (s[i] == 'A')
			d = 'T';
		else if (s[i] == 'T')
			d = 'A';
		else if (s[i] == 'G')
			d = 'C';
		else if (s[i] == 'C') 
			d = 'G';
		
		else if (s[i] == 'a')
			d = 't';
		else if (s[i] == 't')
			d = 'a';
		else if (s[i] == 'g')
			d = 'c';
		else if (s[i] == 'c') 
			d = 'g';
		
		else if (s[i] == '[')
			d = ']';
		else if (s[i] == ']')
			d = '[';
		else 
			d = s[i];
		
		strncat(c, &d, 1);
	}
	
	return c;
	
}

void chomp(char* s) 
{
	
	int j;
	
	j = strlen(s);
	while (j >= 0) {
		if (s[j] == '\n') {
			s[j] = '\0';
			return;
		}
		j--;
	}
	
	
}

char* nextSequence(FILE* fp, char** name, int* size) 
{
	int len = 20000;
	char* seq;
	
	// not yet started, get the first line
	if (nextSequence_started == 0) {
		fgets(nextSequence_currentLine, len, fp);  
		chomp(nextSequence_currentLine);      
		nextSequence_started = 1;
	}
	
	// if file is finished, return 0
	if (nextSequence_ended == 1) {
		return 0;
	}
	
	// if started, line should be filled with ">..."
	if (nextSequence_currentLine[0] == '>') {
		
		*name = strdup(nextSequence_currentLine + 1);
		
		// create a new line
		seq = (char*)calloc(100000, sizeof(char) ); 
		
		while (1) {
			
			fgets(nextSequence_currentLine, len, fp);    
			
			if (feof(fp)) {
				nextSequence_ended  = 1;
				return seq;
			}
			
			chomp(nextSequence_currentLine);      
			if (strlen(nextSequence_currentLine) == 0) { 
				continue;
			}
			
			if (nextSequence_currentLine[0] == '>') {	
				return seq;	
			} else {
				strcat(seq, nextSequence_currentLine);
			}
			
		}
		
	} else return 0;
	
}


char* get_parameter(int argc, char** argv, char* param)
{
	int i = 0;
	while ((i < argc) && (strcmp(param, argv[i])))
		i++;
	if (i<argc)
		return (argv[i+1]);
	else
		return "";
}

//
//  substring
//
char *substr(char *string, int start, int length) 
{
	char *substring;
	
	substring = (char*)calloc(length+1, sizeof(char));
	
	memcpy(substring, string+start, length*sizeof(char));
	substring[length] = '\0';
	
	return substring;
}

