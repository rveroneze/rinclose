#pragma once

#include "stdafx.h"
#include "globalsv.h"
#include "BicsUtils.h"

float runRInCloseCVCP(const dataset_t &D, const row_t &n, const col_t &m, const row_t &minRow, const col_t &minCol, const long &minAttrWeight, const long *attrWeights);
void RInCloseCVCP(const dataset_t &D, const col_t &m, const row_t &minRow, const col_t &minCol, const pbic_t &bic, const long &minAttrWeight, const long *attrWeights, const long *attrWeightsAcc);
bool RCVCP_IsCanonical(const dataset_t &D, const col_t &y, const row_t &p1, const row_t &p2, const pbic_t &bic, col_t &fcol);
