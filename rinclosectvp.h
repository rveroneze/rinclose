#pragma once

#include "stdafx.h"
#include "globalsv.h"
#include "BicsUtils.h"

float runRInCloseCTVP(const dataset_t &D, const row_t &n, const col_t &m, const row_t &minRow, const col_t &minCol);
void RInCloseCTVP(const dataset_t &D, const col_t &m, const row_t &minRow, const col_t &minCol, const pbic_t &bic);
bool RCTVP_IsCanonical(const dataset_t &D, const col_t &y, const row_t &sizeRW, const pbic_t &bic, col_t &fcol);
