#include <stdio.h>
#include <stdlib.h>

#ifdef	__APPLE__
#define	OS	"darwin"
#elif	__linux__
#define OS	"linux"
#else
#define	OS	"unknown"
#endif

/*
 * get the OS
 *
 * currently MacOS X (darwin) and GNU/Linux (linux)
 * are only supported
 */
int
main(int argc, char **argv)
{
	/* print the OS */
	(void)fprintf(stdout, "%s\n", OS);
	
	/* return with success */
	return EXIT_SUCCESS;
}
