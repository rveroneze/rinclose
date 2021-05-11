#pragma once

#include "stdafx.h"

extern unsigned long MVS; // value that represents a Missing Value
extern unsigned short g_output; // 1 - matlab; 2 - python

// ----- Variables for the search of the biclusters -----
extern int g_rnew;
extern unsigned g_cont;
extern pair<data_t, row_t> *g_RWp;
extern ofstream g_filebics;
// -------------------------------------------------------

// ---- Variables for the clique search ----
extern bool **g_Adj;
extern queue<clique_t> g_cliques;
// -----------------------------------------

// ----- Variables for the search of CVC biclusters in OPSM -----
extern row_t *g_RW[2];
// Obs:g_RW[0] is also used RIn-Close_CVC and CTVP
// --------------------------------------------------------------

// ----- Variables for the search using class labels -----
extern unsigned short *g_classes;
extern unsigned short g_maxLabel;
extern double g_minConf;
extern unsigned short g_ignoreLabel;
// --------------------------------------------------------------


// ----- Variable for the experiments about memory RAM usage -----
extern long int g_ram;
// ---------------------------------------------------------------


// ----- Variables for the function RCVC_IsCanonical2 ------------
extern row_t *g_epai;
extern bool *g_rms;
extern data_t *g_menorF, *g_maiorF, *g_menorv, *g_maiorv;
// ---------------------------------------------------------------