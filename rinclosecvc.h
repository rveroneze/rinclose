#pragma once

#include "stdafx.h"
#include "globalsv.h"
#include "BicsUtils.h"

void RCVC_computeRM(const data_t &epsilon, const pbic_t &bic, const pbic_t &child, const row_t &p1, const row_t &p2, const row_t qnmv);

float runRInCloseCVCve(const dataset_t &D, const row_t &n, const col_t &m, const col_t &minCol, const data_t *epsilons);
void RInCloseCVCve(const dataset_t &D, const row_t &n, const col_t &m, const col_t &minCol, const data_t *epsilons, const pbic_t &bic);
bool RCVC_IsCanonical(const dataset_t &D, const data_t *epsilons, const col_t &y, const row_t &p1, const row_t &p2, const pbic_t &bic, col_t &fcol);
bool RCVC_IsCanonical2(const dataset_t &D, const data_t *epsilons, const col_t &y, const row_t &sRW, const pbic_t &bic);