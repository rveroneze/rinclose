#include "globalsv.h"

unsigned long MVS = 999999; // value that represents a Missing Value
unsigned short g_output = 2; // 1 - matlab; 2 - python

// ----- Variables for the search of the biclusters -----
unsigned g_cont = 0; // number of biclusters in the output
pair<data_t, row_t> *g_RWp; // vector to store the pair (data, rows)
ofstream g_filebics; // pointer to the output file
// -------------------------------------------------------

// ----- Variables for the search of CVC biclusters in OPSM -----
row_t *g_RW[2];  // vetor para guardar RW
// Obs:g_RW[0] is also used RIn-Close_CVC
// --------------------------------------------------------------

// ----- Variables for the search using class labels -----
unsigned short *g_classes; // vector to store the class label of each object
unsigned short g_maxLabel; // maximum label
double g_minConf = 0; // confidence threshold
row_t *g_minsups; // vector to store the minsup of each class label
row_t g_smallerMinsup, g_biggerMinsup; // smaller and bigger minsup in *g_minsups
// --------------------------------------------------------------


// ----- Variable for the experiments about memory RAM usage -----
long int g_ram = 0;
// ---------------------------------------------------------------


// ----- Variables for the function RCVC_IsCanonical2 ------------
row_t *g_epai;
bool *g_rms;
data_t *g_menorF, *g_maiorF, *g_menorv, *g_maiorv;
// ---------------------------------------------------------------