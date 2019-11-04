#include "stdafx.h"
#include "globalsv.h"
#include "IK_GPX.h"
#include "BicsUtils.h"

void S3_step3(const dataset_t &Da, const row_t &n, const col_t &mAug, const col_t &minCol2, const data_t &epsilon, const pbic_t &bic);
bool S3_IsMaximal(const dataset_t &Da, const row_t &n, const col_t &m, const data_t &epsilon, const pbic_t &bic, const clique_t &clique, const col_t *B2);
col_t S3_getColDa(const col_t &m, const col_t &col1, const col_t &col2);
void S3_getColsOrig(const col_t &m, const col_t &colAug, col_t &c1, col_t &c2);

void S3_step3_OP(const datasetOP_t &Da, const row_t &n, const col_t &mAug, const col_t &minCol2, const pbic_t &bic);
bool S3_IsMaximal_OP(const datasetOP_t &Da, const row_t &n, const col_t &m, const pbic_t &bic, const clique_t &clique, const col_t *B2);