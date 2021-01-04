#pragma once

#include "stdafx.h"
#include "globalsv.h"

void openPrintFile(const string &filename);
void printBic(const pbic_t &bic, const col_t m);
void closePrintFile();
