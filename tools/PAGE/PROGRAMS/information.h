/* information.c */

typedef struct _FloatAndIndex {
  int   i;
  float v;
} FloatAndIndex;
  


typedef struct _MI_Contribs {
  
  int i;
  int j;
  float co;
  
} MI_Contribs;

int isDiagonalDominant(int* A_q, int* B_q, int DA, int DB, int n); 

int max_rank_test_cond(double score, int* M_q, int mbins, int* E_q, int ebins, int* A_q, int abins, int nborfs, int shuffle, double* max_cmi_shu, int* val, double* z); 
double InformationIncrease(int* A_q, int DA, int* B_q, int DB, int* E_q, int DE, int n, double* I1, double* I2); 

int CmpFI(const void* _a, const void* _b);

void StudyContribution_FindInfoLocalMI(int** counts, int D1, int D2, MI_Contribs** mico); 
void InformationGain_D2(int** counts, int D1, int D2, MI_Contribs** mico, int sort); 

float CalculateCondMIbasic( int* data1,  int* data2, int* data3, int m, int D1, int D2, int D3); 
int** Construct3DCountMatrix(int* qv1, int* qv2, int* qv3, int m, int D1, int D2, int D3);
float CalculateNormalizedMIbasic(int* data1,  int* data2, int m, int D1, int D2); 


double FindCondMI(int*** counts, int D1, int D2, int D3); 

float  CalculateMIbasic(int* data1,  int* data2, int m, int D1, int D2); 
void   QuantizePointerized(float* v, int n, int D, int** qv, float** bins); 

float* sortFloatVector(float* v, int n);
int*   ThresholdQuantizeInt(int* v, int n, int t); 

int*   TopQuantize(float* v, int n, int ne, float* th);
int*   ThresholdQuantize(float* v, int n, float t, int sup); 

int**  allocateQuantizedMatrix(int n); 

float  CalculateMIFromUniquePreQuantizedVectors(int* quantized_vector1, int D1, int* quantized_vector2, int D2, int m, int Repeats, int remove_missing_values, int correct_sample_size); 

float entropy_from_probas(float* P, int n); 

float FindInfo(int *qv1, int *qv2, int FullSample, int D1, int D2, float *SampleSizesPerc, int N, int *Trials, int Repeats, float* std);
int **ConstructCountMatrix(int *qv1, int *qv2, int m, int D1, int D2);
void freeCountMatrix(int **counts, int D1, int D2);
double FindInfoLocalMI(int **counts, int D1, int D2);

int* Quantize(float* v, int n, int D, float** bins); 

int CmpDbl(const void *_a, const void *_b);
void showCountMatrix(int** counts, int D1, int D2, char** states1, char** states2, char* prefix); 

float QuantizeAndCalculateMI(float* data1, float* data2, int m, int* quantized_data1, int QD1, int* quantized_data2, int QD2); 

float CalculateMI(int* data1,  int* data2, int m, int D1, int D2, int Repeats); 
int* combineQuantizedVectors(int* data1,  int* data2, int m, int D1, int D2); 
int* combineBinaryVectorsAND(int* data1,  int* data2, int m); 

void CalculateAndShowCountMatrix(int* data1,  int* data2, int m, int D1, int D2, char** states1, char** states2, char* prefix); 

void  removeMissingValuesFloatFloat(float* data1, float* data2, 
				    int m, 
				    float** data1_nomissing, 
				    float** data2_nomissing,
				    int* m_nomissing); 
void  removeMissingValuesIntFloat(int* data1, float* data2, 
				 int m, 
				 int** data1_nomissing, 
				 float** data2_nomissing,
				 int* m_nomissing); 
void  removeMissingValuesIntInt(int* data1, int* data2, 
				int m, 
				int** data1_nomissing, 
				int** data2_nomissing,
				int* m_nomissing); 
float CalculateMIFromPreQuantizedVectors(int** quantized_vectors1, int* vqd1, int nqd1, int** quantized_vectors2, int* vqd2, int nqd2, int m, int* OD1, int* OD2, int Repeats, int remove_missing_values, int correct_sample_size); 
void  removeMissingValuesFloat(float* data1, int m, float** data1_nomissing, int* m_nomissing, char** mask); 
void  getMaskedVector(int* data1, int m, char* mask, int** data1_nomissing); 


float KL_dist(float* P1, float* P2, int D);
float JS_dist(float* P1, float* P2, int D); 
 

float entropy(int* data, int m, int D); 
int* logicCombine(int* data1,  int* data2, int m, int D1, int D2, char* type); 
