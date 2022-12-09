#include "stdafx.h"
#include "globalsv.h"

unsigned long MVS = 999999; // value that represents a Missing Value
unsigned short g_output = 2; // 1 - matlab; 2 - python

// ----- Variables for the search of the biclusters -----
unsigned g_cont = 0; // number of biclusters in the output
pair<data_t, row_t> *g_RWp; // vector to store the pair (data, rows)
ofstream g_filebics; // pointer to the output file
// -------------------------------------------------------

// ---- Variables for the clique search - CHV biclusters ----
bool **g_Adj; // stores the adjacency matrix
queue<clique_t> g_cliques; // stores cliques
// -----------------------------------------

// ----- Variables for the search of CVC biclusters in OPSM -----
row_t *g_RW[2];  // vetor para guardar RW
// Obs:g_RW[0] is also used RIn-Close_CVC and CTVP
// --------------------------------------------------------------

// ----- Variables for the search using class labels -----
unsigned short *g_classes; // vector to store the class label of each object
unsigned short g_maxLabel; // maximum label
double g_minConf = 0; // confidence threshold
unsigned short g_ignoreLabel; // biclusters that represents label g_ignoreLabel will not be outputted
// --------------------------------------------------------------


// ----- Variable for the experiments about memory RAM usage -----
long int g_ram = 0;
// ---------------------------------------------------------------


// ----- Variables for the function RCVC_IsCanonical2 ------------
row_t *g_epai;
bool *g_rms;
data_t *g_menorF, *g_maiorF, *g_menorv, *g_maiorv;
// ---------------------------------------------------------------