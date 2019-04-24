#include <errno.h>
#define __set_errno(x) errno=x

#include <sys/types.h>

#include <string.h>
#include <search.h>
#include <stdlib.h>
#include <stdio.h>

#include "hashtable.h"

  
void HASH_init(HASH* h, int n)
{
  int hashret;
  h->hash_table = (struct my_hsearch_data*)calloc(1,  sizeof(struct my_hsearch_data));
  hashret = my_hcreate_r(n, h->hash_table);
  if (hashret == 0) {
    printf("Could not create hash table with %d elements ...\n", n);
    exit(0);
  } 
}

void HASH_insert(HASH* h, char* k, int v)
{
  HASH_enter(h, k, v);
}

void HASH_enter(HASH* h, char* k, int v)
{
  int hashret;
  (h->e).key   = strdup( k );
  (h->e).data  = (char*)v;
  hashret = my_hsearch_r(h->e, ENTER, &(h->ep), h->hash_table);
  if (hashret == 0) {
    printf("Could not enter entry into hash table ...\n");
    exit(0);
  }
}


int HASH_find(HASH* h, char* k, int* v)
{
  int hashret;
  (h->e).key   = k;
  hashret = my_hsearch_r(h->e, FIND, &(h->ep), h->hash_table);
  //if (hashret == 0) {
  //  printf("Could not FIND entry into hash table ...\n");
  //  exit(0);
  //}
  if (h->ep) {
    *v = (int)(h->ep->data);
    return 1;
  } else {
    *v = 0;
    return 0;
  }
}


void HASH_destroy(HASH* h)
{
  my_hdestroy_r(h->hash_table);
}




/* [Aho,Sethi,Ullman] Compilers: Principles, Techniques and Tools, 1986
   [Knuth]            The Art of Computer Programming, part 3 (6.4)  */


/* The reentrant version has no static variables to maintain the state.
   Instead the interface of all functions is extended to take an argument
   which describes the current status.  */

/* For the used double hash method the table size has to be a prime. To
   correct the user given table size we need a prime test.  This trivial
   algorithm is adequate because
   a)  the code is (most probably) called a few times per program run and
   b)  the number is small because the table must fit in the core  */

static int my_isprime (unsigned int number)
{
  /* no even number will be passed */
  unsigned int div = 3;

  while (div * div < number && number % div != 0)
    div += 2;

  return number % div != 0;
}


/* Before using the hash table we must allocate memory for it.
   Test for an existing table are done. We allocate one element
   more as the found prime number says. This is done for more effective
   indexing as explained in the comment for the hsearch function.
   The contents of the table is zeroed, especially the field used
   becomes zero.  */
int
my_hcreate_r (size_t nel, struct my_hsearch_data *htab)
{
  /* Test for correct arguments.  */
  if (htab == NULL)
    {
      __set_errno (EINVAL);
      return 0;
    }

  /* There is still another table active. Return with error. */
  if (htab->table != NULL)
    return 0;

  /* Change nel to the first prime number not smaller as nel. */
  nel |= 1;      /* make odd */
  //printf("About to enter isprime\n"); fflush(stdout);
  while (!my_isprime (nel)) {
    //printf("trying nel=%d\n", nel); fflush(stdout);
    nel += 2;
  }

  htab->size = nel;
  htab->filled = 0;

  /* allocate memory and zero out */
  htab->table = (_ENTRY *) calloc (htab->size + 1, sizeof (_ENTRY));
  if (htab->table == NULL)
    return 0;

  /* everything went alright */
  return 1;
}



/* After using the hash table it has to be destroyed. The used memory can
   be freed and the local static variable can be marked as not used.  */
void
my_hdestroy_r (struct my_hsearch_data *htab)
{
  /* Test for correct arguments.  */
  if (htab == NULL)
    {
      __set_errno (EINVAL);
      return;
    }

  if (htab->table != NULL)
    /* free used memory */
    free (htab->table);

  /* the sign for an existing table is an value != NULL in htable */
  htab->table = NULL;
}




/* This is the search function. It uses double hashing with open addressing.
   The argument item.key has to be a pointer to an zero terminated, most
   probably strings of chars. The function for generating a number of the
   strings is simple but fast. It can be replaced by a more complex function
   like ajw (see [Aho,Sethi,Ullman]) if the needs are shown.

   We use an trick to speed up the lookup. The table is created by hcreate
   with one more element available. This enables us to use the index zero
   special. This index will never be used because we store the first hash
   index in the field used where zero means not used. Every other value
   means used. The used field can be used as a first fast comparison for
   equality of the stored and the parameter value. This helps to prevent
   unnecessary expensive calls of strcmp.  */
int
my_hsearch_r (
	      ENTRY item,
	      ACTION action,
	      ENTRY **retval,
	      struct my_hsearch_data *htab)
{
  unsigned int hval;
  unsigned int count;
  unsigned int len = strlen (item.key);
  unsigned int idx;

  /* Compute an value for the given string. Perhaps use a better method. */
  hval = len;
  count = len;
  while (count-- > 0)
    {
      hval <<= 4;
      hval += item.key[count];
    }

  /* First hash function: simply take the modul but prevent zero. */
  hval %= htab->size;
  if (hval == 0)
    ++hval;

  /* The first index tried. */
  idx = hval;

  if (htab->table[idx].used)
    {
      /* Further action might be required according to the action value. */
      unsigned hval2;

      if (htab->table[idx].used == hval
	  && strcmp (item.key, htab->table[idx].entry.key) == 0)
	{
	  *retval = &htab->table[idx].entry;
	  return 1;
	}

      /* Second hash function, as suggested in [Knuth] */
      hval2 = 1 + hval % (htab->size - 2);

      do
	{
	  /* Because SIZE is prime this guarantees to step through all
             available indeces.  */
          if (idx <= hval2)
	    idx = htab->size + idx - hval2;
	  else
	    idx -= hval2;

	  /* If we visited all entries leave the loop unsuccessfully.  */
	  if (idx == hval)
	    break;

            /* If entry is found use it. */
          if (htab->table[idx].used == hval
	      && strcmp (item.key, htab->table[idx].entry.key) == 0)
	    {
	      *retval = &htab->table[idx].entry;
	      return 1;
	    }
	}
      while (htab->table[idx].used);
    }

  /* An empty bucket has been found. */
  if (action == ENTER)
    {
      /* If table is full and another entry should be entered return
	 with error.  */
      if (htab->filled == htab->size)
	{
	  __set_errno (ENOMEM);
	  *retval = NULL;
	  return 0;
	}

      htab->table[idx].used  = hval;
      htab->table[idx].entry = item;

      ++htab->filled;

      *retval = &htab->table[idx].entry;
      return 1;
    }

  __set_errno (ESRCH);
  *retval = NULL;
  return 0;
}

