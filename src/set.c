#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "set.h"

/* item comparator; internal */
static int
cmp(const void *lhs, const void *rhs)
{
	/* wrap strcmp() */
	return strcmp(*(char **)lhs, *(char **)rhs);
}

/* 
 * initialize the set
 *
 * @s:		the set to initialize
 * returns:	0 on success, 1 on failure
 */
int
set_init(set_t *s)
{
	/* sanity checks */
	if (s == NULL)
		/* failed */
		return 1;

	/* allocate initial set space */
	if ((s->arr = calloc(ALLOC_BLOCK, sizeof(char *))) == NULL)
		/* failed */
		return 1;

	/* set the capacity */
	s->capacity = ALLOC_BLOCK;

	/* set the size */
	s->size = 0;
	
	/* success */
	return 0;
}

/*
 * add item to the set
 *
 * @s:		the set to add the item
 * @item:	the item to be added
 * returns:	0 on success, -1 on failure
 */
int
set_add_item(set_t *s, char *item)
{
	char *fit 	= NULL;	/* found flag; initially NULL */
	char **narr	= NULL;	/* new array; for realloc() */
	/* sanity checks */
	if (s == NULL || item == NULL)
		/* failed */
		return 1;

	/* set is full; time to expand */
	if (s->size == s->capacity) {

		/* try to reallocate */
		if ((narr = realloc(s->arr,
			(s->capacity + ALLOC_BLOCK) * sizeof(char *))) == NULL)
			/* failed */
			return 1;

		/* patch the pointers */
		s->arr = narr;
		
		/* update the capacity */
		s->capacity += ALLOC_BLOCK;
	}

	/* 
	 * check to see if the item already exists on the set;
	 * it uses binary search, therefore this should be O(N lg N) 
	 */
	if ((fit = bsearch(&item,
		s->arr,
		s->size,
		sizeof(char *),
		cmp)) == NULL) {
			
		/* not found; add it */
		s->arr[s->size] = item;
		s->size++;
			
		/* sort the set; this should be N lg N on average */
		qsort(s->arr, s->size, sizeof(char *), cmp);

		/* just added */
		return 0;
	}
	else
		/* already in the set */
		return 1;
}

/*
 * dispose a set
 *
 * s:	the set to dispose
 */
void
set_del(set_t *s)
{
	/* clearup */
	if (s->arr != NULL)
		free(s->arr);
}
