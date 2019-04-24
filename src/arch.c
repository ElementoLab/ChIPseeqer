#include <stdio.h>
#include <stdlib.h>

#ifdef	__i386__
#define	ARCH	"i386"
#elif	__x86_64__
#define ARCH	"x86_64"
#else
#define	ARCH	"unknown"
#endif

/*
 * get the system architecture
 *
 * currently i386 and x86_64 are
 * only supported
 */
int
main(int argc, char **argv)
{
	/* print the system arch */
	(void)fprintf(stdout, "%s\n", ARCH);
	
	/* return with success */
	return EXIT_SUCCESS;
}
