#pragma once

#include "stdafx.h"
#include "globalsv.h"
//#include "BicsUtils.h"
#include "step3.h"

float runRInCloseCVCP_OP(const datasetOP_t &D, const row_t &n, const col_t &m, const row_t &minRow, const col_t &minCol);
void RInCloseCVCP_OP(const datasetOP_t &D, const row_t &n, const col_t &m, const row_t &minRow, const col_t &minCol, const pbic_t &bic);
bool RCVCP_OP_IsCanonical(const datasetOP_t &D, const col_t &y, const row_t &p1, const row_t &p2, const pbic_t &bic, col_t &fcol);
