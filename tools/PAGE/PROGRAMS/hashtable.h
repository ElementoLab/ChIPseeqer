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

int  my_hcreate_r  (size_t nel,  struct my_hsearch_data *htab);
void my_hdestroy_r (struct my_hsearch_data *htab);
int  my_hsearch_r  (ENTRY item, ACTION action, ENTRY **retval, struct my_hsearch_data *htab);
     

