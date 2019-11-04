#pragma once

#include "stdafx.h"
#include "globalsv.h"

void openPrintFile(const string &filename);
void printBic(const pbic_t &bic, const col_t m);
void printBic(const pbic_t &bic, const col_t *B, const col_t sizeB);
void printBic(const pbic_t &bic, const clique_t &clique, const col_t *B);
void closePrintFile();
double getMinConf(row_t *A, row_t size);
