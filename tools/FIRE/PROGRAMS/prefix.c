
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prefix.h"


void PT_allocateMoreNode(PT* t, int n) 
{

  int t_size_tmp;
  prefixNode* t_pnode_tmp;
  int* t_node2index_tmp;

  t_size_tmp = t->size + n;
  
  t_node2index_tmp = (int*)malloc(t_size_tmp * sizeof(int));
  memcpy(t_node2index_tmp, t->node2index, t->size * sizeof(int));
  free(t->node2index);
  t->node2index    = t_node2index_tmp;
  
  t_pnode_tmp = (prefixNode*)malloc(t_size_tmp * (sizeof(prefixNode)));
  memcpy(t_pnode_tmp, t->pnode, t->size * sizeof(prefixNode));  
  free(t->pnode);
  t->pnode    = t_pnode_tmp;

  t->size = t_size_tmp;
}


void PT_createPrefixTree(PT* t, int n, int motifsize) 
{
  
  t->size = n * motifsize;
  
   // associates a node to an index
  t->node2index = (int*)malloc(t->size * sizeof(int));
  if (t->node2index == 0) {
    printf("Not enough memory for prefix tree (node2index)\n");
  }


  // we are going to need many nodes ...
  t->pnode = (prefixNode*)malloc(t->size * (sizeof(prefixNode)));
  if (t->pnode == 0) {
    printf("Not enough memory for prefix tree (pnode)\n");
  }

  (t->pnode)[0].next[0] = -1;
  (t->pnode)[0].next[1] = -1;
  (t->pnode)[0].next[2] = -1;
  (t->pnode)[0].next[3] = -1;
  //(t->pnode)[0].next[4] = -1;


  //
  //  add all the words in the list
  //
  t->nextfree  = 1;
  t->nodes_to_add = 100000;
  t->n_allowed = 0;
  
}


void PT_populate(PT* pt, char** kmers, int nbkmers) {
  int i;
  int lastnode;

  for (i=0; i<nbkmers; i++) {
    lastnode  =  PT_addWord(pt, kmers[i]);
    (pt->node2index)[lastnode] = i;
  }
  
}


// add word to the prefix tree (if it is not there)
int PT_existWord(PT* t, char* s) 
{

  int i, k;
  int tt;
  int j = 0;
  int m = strlen(s);
  int cnt_n = 0;
  
  for (i=0; i<m; i++) {
    if (s[i] == 'N')
      cnt_n ++;
  }
  if (cnt_n > 2 * t->n_allowed)
    return -2;
  
  for (i=0; i<m; i++) {
    tt = nucl(s[i]);
    
    if (tt == -1)
      return -2;    // case where character is not allowed

    if ((s[i] == 'N') && (t->n_allowed == 0))
      return -2;

    if ((s[i] == 'N') && (i+1 > t->n_allowed) && (i < m - t->n_allowed))
      return -2;  

    k = (t->pnode)[j].next[nucl(s[i])];
    if (k == -1)
      return -1;
  
    j = (t->pnode)[j].next[nucl(s[i])];

  }

  return j;
}



// add word to the prefix tree (if it is not there)
int PT_addWord(PT* t, char* s) 
{
  int i, k;
  int j = 0;
  int m = strlen(s);
  for (i=0; i<m; i++) {
    k = (t->pnode)[j].next[nucl(s[i])];
    if (k == -1) {
      (t->pnode)[j].next[nucl(s[i])] = t->nextfree;
      (t->pnode)[t->nextfree].next[0] = -1;
      (t->pnode)[t->nextfree].next[1] = -1;
      (t->pnode)[t->nextfree].next[2] = -1;
      (t->pnode)[t->nextfree].next[3] = -1;
      //(t->pnode)[t->nextfree].next[4] = -1;
      (t->nextfree)++;

      if (t->nextfree == t->size) {
	printf("Warning: prefix tree full, allocating more nodes.\n");
	PT_allocateMoreNode(t, t->nodes_to_add); 
      }

    }
    j = (t->pnode)[j].next[nucl(s[i])];
  }
  t->lastnode = j;
    
  return j;
}







void InitAndCreatePrefixTree(prefixNode** pnode, int** node2index, char*** kmers, int n, int motifsize) 
{

  int nextfree;
  int i;
  int lastnode;
  
  int size = n * motifsize;
  
   // associates a node to an index
  *node2index = (int*)malloc(size * sizeof(int));
  
  // we are going to need many nodes ...
  *pnode = (prefixNode*)malloc(size * (sizeof(prefixNode)));
  (*pnode)[0].next[0] = -1;
  (*pnode)[0].next[1] = -1;
  (*pnode)[0].next[2] = -1;
  (*pnode)[0].next[3] = -1;


  //
  //  add all the words in the list
  //
  nextfree = 1;

  for (i=0; i<n; i++) {
    lastnode = addWord((*kmers)[i], motifsize, *pnode, &nextfree);

    //
    //   at this point, lastnode contains the last node inserted in the tree
    //
    
    (*node2index)[lastnode] = i;
    
    
    
  }
  


}


// add word to the prefix tree (if it is not there)
int existWord(char* s, int m, prefixNode* pnode) 
{

  //printf("Adding %s to the prefix tree ..\n", s);
  int i, k;
  int t;
  int j = 0;

  for (i=0; i<m; i++) {
    t = nucl(s[i]);
    
    if (t == -1)
      return -2;

    k = pnode[j].next[nucl(s[i])];
    if (k == -1)
      return -1;
  
    j = pnode[j].next[nucl(s[i])];

  }

  return j;
}



// add word to the prefix tree (if it is not there)
int addWord(char* s, int m, prefixNode* pnode, int* nextfree) 
{

  //printf("Adding %s to the prefix tree ..\n", s);
  int i, k;

  int j = 0;

  for (i=0; i<m; i++) {
    k = pnode[j].next[nucl(s[i])];
    if (k == -1) {
      pnode[j].next[nucl(s[i])] = *nextfree;

      pnode[*nextfree].next[0] = -1;
      pnode[*nextfree].next[1] = -1;
      pnode[*nextfree].next[3] = -1;
      pnode[*nextfree].next[2] = -1;

      (*nextfree)++;
    }
    j = pnode[j].next[nucl(s[i])];
  }

  return j;
}



int nucl(char c) 
{
  
  if (c == 'A')
    return 0;
  else if (c == 'C')
    return 1;
  else if (c == 'G')
    return 2;
  else if (c == 'T')
    return 3;
  //  else if (c == 'N')
  //  return 4;
  else return -1;
  
  
}



