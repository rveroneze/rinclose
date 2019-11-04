#pragma once

#include "stdafx.h"
//#include "rinclosecvc_op.h"  // versao usando o cvc - nem implementei - se precisar: pegar da versao anterior
#include "rinclosecvcp_op.h"  // versao usando o cvcp
#include "step3.h"

float runRInCloseOPSM(const dataset_t &D, const row_t &n, const col_t &m, const row_t &minRow, const col_t &minCol);
