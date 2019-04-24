

typedef struct _prefixNode {
  int next[4];    // A, C, G, T, N
 } prefixNode;

typedef struct _PT {
  prefixNode*  pnode;
  int*         node2index;
  int          nextfree;
  int          lastnode;
  int          size;
  int          n_allowed;
  int          nodes_to_add;
} PT;



void PT_createPrefixTree(PT* t,  int n, int motifsize); 
int  PT_existWord       (PT* t,  char* s); 
int  PT_addWord         (PT* t,  char* s); 
void PT_populate        (PT* pt, char** kmers, int nbkmers);


int nucl(char c); 

int addWord(char* s, int m, prefixNode* pnode, int* nextfree); 
int existWord(char* s, int m, prefixNode* pnode); 
void InitAndCreatePrefixTree(prefixNode** pnode, int** node2index, char*** kmers, int n, int motifsize); 

