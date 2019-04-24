#include <errno.h>
#include <getopt.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define BUF_SZ		128			/* buffer size */

#define VERSION	"1.0"				/* version number */
#define AUTHOR 	"Eugenia G. Giannopoulou"	/* author info */
#define EMAIL	"eug2002@med.cornell.edu"	/* author email */

/* getopt companions */
extern char *optarg;
extern int optind, opterr, optopt;

/* function prototypes */
static void fill_indx(size_t**, size_t, size_t, char**);
static int cmp(const void*, const void*);
static void help(char **argv);
static void version(char **version);

/* time 0, place 0 */
int
main(int argc, char **argv) {
	
	/* getopt stuff */
	char opt;				/* option */
	int long_opt_index		= 0;    /* indx */
	char *input			= NULL;	/* input filename */
	int ln				= 0;	/* input file size
						   (total lines number) */
	FILE *f				= NULL;	/* input file */
	int rnd				= 0;	/* random lines number */
	size_t *indx			= NULL;	/* rnd numbers array */
	char line_buf[BUF_SZ];			/* line buffer */
	size_t i			= 0;	/* iterator */
	size_t indx_i			= 0;	/* rnd numbers array iterator */

	/* the long options */
	struct option long_options[] = {
		{"file",        1, NULL, 'f'},  /* -f / --file */
		{"size",	1, NULL, 's'},	/* -s / --size */
		{"lines",       1, NULL, 'l'},  /* -l / --lines */
		{"help",        0, NULL, 'h'},  /* -h / --help */
		{"version",     0, NULL, 'v'},  /* -v / --version */
		{0, 0, 0, 0}                    /* terminating item */
	};

	opterr			= 0;    /* no default error messages */

	/* parse input */
	while ((opt = getopt_long(argc, argv, ":hvf:l:s:", long_options,
						&long_opt_index)) != -1) {
                switch (opt) {
			/* input file */
                        case 'f':
				input = optarg;
				break;
			/* input file size */
                        case 's':
				ln = atoi(optarg);
				break;
			/* linux number */
			case 'l':
				rnd = atof(optarg);
				break;
			/* help */
                        case 'h':
				help(argv);
			/* version */
                        case 'v':
				version(argv);
				exit(0);
				break;
			/* unknown option */
			case '?':
				(void)fprintf(stderr, "%s: invalid option -- %c\n",
								argv[0], (char)optopt);
				(void)fprintf(stderr, "Try `%s --help for ", argv[0]);
				(void)fprintf(stderr, "more information\n");
				exit(1);
				break;  /* not reached */
			/* missing argument */
			case ':':
				(void)fprintf(stderr, "%s: option requires an argument -- %c\n",
								argv[0], (char)optopt);
				(void)fprintf(stderr, "Try `%s --help for ", argv[0]);
				(void)fprintf(stderr, "more information\n");
				exit(2);
                                break;  /* not reached */
			/* make the compiler happy */
                        default:        /* not reached */
                                break;  /* not reached */

		}
	}

	/* input validation */
	if (ln <= 0) {
		/* failed */
		(void)fprintf(stderr, "%s: invalid total number of lines argument -- %d\n", argv[0], ln);
		(void)fprintf(stderr, "Try `%s --help for ", argv[0]);
		(void)fprintf(stderr, "more information\n");
		exit(3);
	}
	if (input == NULL) {
		/* failed */
		(void)fprintf(stderr, "%s: invalid file argument -- %s\n", argv[0], input);
		(void)fprintf(stderr, "Try `%s --help for ", argv[0]);
		(void)fprintf(stderr, "more information\n");
		exit(4);
	}
	if ((rnd <= 0) | (rnd > ln)) {
		/* failed */
		(void)fprintf(stderr, "%s: invalid lines argument -- %d\n", argv[0], rnd);
		(void)fprintf(stderr, "Try `%s --help for ", argv[0]);
		(void)fprintf(stderr, "more information\n");
		exit(5);
	}
	
	/* open the input file */
	if ((f = fopen(input, "r")) == NULL) {
		/* failed */
		(void)fprintf(stderr, "%s: failed while openning input file -- %s\n", argv[0], 
										strerror(errno));
		(void)fprintf(stderr, "Try `%s --help for ", argv[0]);
		(void)fprintf(stderr, "more information\n");
		exit(6);
	}
	
	/* fll the rnd array */
	fill_indx(&indx, rnd, (size_t)ln, argv);
	
	/* read the file line-by-line */
	for(i=1; i<=ln; i++) {

		/* clean the line buffer */
		(void)memset(line_buf, 0, BUF_SZ);

		/* read a line */
		(void)fgets(line_buf, BUF_SZ, f);

		/* split */
		if (indx_i < rnd && indx[indx_i] == i) {
			/* print to stdout */
			(void)fprintf(stdout, "%s", line_buf);
			/* advance the index */
			indx_i++;
		}
		/* else */
			/* print to stdout */
			/* (void)fprintf(stdout, "%s", line_buf); */
	}

	/* cleanup */
	(void)fclose(f);
	free(indx);
	
	/* bye bye */
	return EXIT_SUCCESS;
}


/* fill the rnd numbers array */
static void
fill_indx(size_t **r, size_t rln, size_t ln, char **argv) {
	
	size_t i;				/* iterator */
	size_t norm_fact	= 1;		/* normalizer */
	double norm_val 	= 0.0;		/* normalized value */
	size_t rndln		= 0;		/* random line number */

	/* fix normalizer to 2**31 - 1 */
	norm_fact <<= 31;
	norm_fact -= 1;

	/* init the pseudo-random generator */
	srandom((unsigned)(getpid() + time(NULL)));

	/* alloc space */
	if ((*r = malloc(rln * sizeof(size_t))) == NULL) {
		/* failed */
		(void)fprintf(stderr, "%s: failed while allocating memory -- %s\n", argv[0], 
									strerror(errno));
		(void)fprintf(stderr, "Try `%s --help for ", argv[0]);
		(void)fprintf(stderr, "more information\n");
	}

	/* init the array */
	for (i = 0; i < rln; i++)
		(*r)[i] = ln + 1;

	/* fill the array */
	for(i = 0; i < rln; i++) {
again:
		/* normalize the random value to [0..1] */
		norm_val = (double)random() / norm_fact;

		/* get a value to [1..ln] */
		rndln = (size_t)lrintl(((ln - 1) * norm_val) + 1);
	
		/* check if it exists on the array - use binary search */
		if (bsearch(&rndln, *r, rln, sizeof(size_t), cmp) == NULL) 
			/* not found - put at the end */
			(*r)[rln - 1] = rndln;
		else
			/* found - try again */
			goto again;
		
		/* sort the array - quicksort */
		qsort(*r, rln, sizeof(size_t), cmp);
	}	
}

/* compare function - used for searching and sorting */
static int
cmp(const void *a, const void *b) {
	return *((size_t*)a) - *((size_t*)b);
}

/*
 * help info
 */
static void
help(char **argv) {
	(void)fprintf(stderr, "Usage: %s [OPTION]\n", argv[0]);
	(void)fprintf(stderr, "Executes the split utility.\n\n");
	(void)fprintf(stderr, "  -f, --file\tthe file to split\n"); 
	(void)fprintf(stderr, "  -s, --size\ttotal number of lines in file\n"); 
	(void)fprintf(stderr, "  -l, --lines\trandom lines number\n"); 
	(void)fprintf(stderr, "  -h, --help\tdisplay this help "); 
	(void)fprintf(stderr, "and exit\n");
	(void)fprintf(stderr, "  -v, --version\toutput version "); 
	(void)fprintf(stderr, "information and exit\n\n");
	(void)fprintf(stderr, "  Report bugs to <%s>.\n\n", EMAIL);
}

/*
 * version info
 */
static void
version(char **argv) {
	(void)fprintf(stderr, "%s %s\n", argv[0], VERSION);
	(void)fprintf(stderr, "Copyright (C) 2009 by %s\n\n", AUTHOR);
	(void)fprintf(stderr, "Permission to use, copy, and modify this "); 
	(void)fprintf(stderr, "software without fee is hereby\ngranted, ");
	(void)fprintf(stderr, "provided that this entire notice is included "); 
	(void)fprintf(stderr, "in all copies of any\nsoftware which is or ");
	(void)fprintf(stderr, "includes a copy or modification of this software.\n\n");
	(void)fprintf(stderr, "THIS SOFTWARE IS BEING PROVIDED \"AS IS\", ");
	(void)fprintf(stderr, "WITHOUT ANY EXPRESS OR IMPLIED\nWARRANTY. ");
	(void)fprintf(stderr, "IN PARTICULAR, THE AUTHOR MAKES NO ");
	(void)fprintf(stderr, "REPRESENTATION OR WARRANTY OF\nANY KIND "); 
	(void)fprintf(stderr, "CONCERNING THE MERCHANTABILITY OF THIS ");
	(void)fprintf(stderr, "SOFTWARE OR ITS FITNESS\nFOR ANY PARTICULAR ");
	(void)fprintf(stderr, "PURPOSE\n\n");
}
