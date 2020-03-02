// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include <stdio.h>


// TODO: reference additional headers your program requires here
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <time.h> 
#include <algorithm>
#include <math.h>
#include <queue>
#include <float.h>
#include <unordered_map>  // usado para o arquivo config.txt

#include <sys/resource.h> // for the experiments about memory RAM usage

using namespace std;

typedef double data_t;
typedef data_t** dataset_t;

typedef unsigned int dataOP_t;
typedef dataOP_t** datasetOP_t;

typedef unsigned int row_t;
typedef unsigned int col_t;

struct bic_t {
	row_t *A;
	row_t sizeA;
	bool *B, *PN;
	col_t sizeB;
	row_t sizeRM;
	col_t col;
	long AttrWeight;
};

typedef bic_t *pbic_t;

struct clique_t {
	col_t size;
	col_t *nodes;
};
