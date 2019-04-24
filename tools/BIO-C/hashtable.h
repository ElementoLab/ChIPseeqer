#ifndef HASHTABLE_H
#define HASHTABLE_H 

typedef struct _ENTRY
{
  unsigned int used;
  ENTRY entry;
}
_ENTRY;


struct my_hsearch_data
{
  int filled;
  _ENTRY* table;
  int size;
}; 

typedef struct _HASH {
  ENTRY e;
  ENTRY* ep;
  struct my_hsearch_data* hash_table;
} HASH;


int  my_hcreate_r  (size_t nel,  struct my_hsearch_data *htab);
void my_hdestroy_r (struct my_hsearch_data *htab);
int  my_hsearch_r  (ENTRY item, ACTION action, ENTRY **retval, struct my_hsearch_data *htab);
     

void HASH_init(HASH* h, int n);
void HASH_enter(HASH* h, char* k, int v);
void HASH_insert(HASH* h, char* k, int v);
int  HASH_find(HASH* h, char* k, int* v);
void HASH_destroy(HASH* h);


#endif
