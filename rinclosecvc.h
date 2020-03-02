#pragma once

#include "stdafx.h"
#include "globalsv.h"
#include "BicsUtils.h"
#include "step3.h"

float runRInCloseCVC(const dataset_t &D, const row_t &n, const col_t &m, const row_t &minRow, const col_t &minCol, const data_t &epsilon, const bool &s3);
void RInCloseCVC(const dataset_t &D, const row_t &n, const col_t &m, const row_t &minRow, const col_t &minCol, const data_t &epsilon, const pbic_t &bic, const bool &s3);
bool RCVC_IsCanonical(const dataset_t &D, const data_t &epsilon, const col_t &y, const row_t &p1, const row_t &p2, const pbic_t &bic, col_t &fcol);
bool RCVC_IsCanonical2(const dataset_t &D, const data_t &epsilon, const col_t &y, const row_t &sRW, const pbic_t &bic);
void RCVC_computeRM(const row_t &minRow, const data_t &epsilon, const pbic_t &bic, const pbic_t &child, const row_t &p1, const row_t &p2, const row_t qnmv);

float runRInCloseCVCve(const dataset_t &D, const row_t &n, const col_t &m, const row_t &minRow, const col_t &minCol, const data_t *epsilons, const long &minAttrWeight, const long *attrWeights);
void RInCloseCVCve(const dataset_t &D, const row_t &n, const col_t &m, const row_t &minRow, const col_t &minCol, const data_t *epsilons, const pbic_t &bic, const long &minAttrWeight, const long *attrWeights, const long *attrWeightsAcc);
bool RCVC_IsCanonical(const dataset_t &D, const data_t *epsilons, const col_t &y, const row_t &p1, const row_t &p2, const pbic_t &bic, col_t &fcol);
bool RCVC_IsCanonical2(const dataset_t &D, const data_t *epsilons, const col_t &y, const row_t &sRW, const pbic_t &bic);