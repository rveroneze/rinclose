#pragma once

#include "stdafx.h"
#include "globalsv.h"

void runIK_GPX(const col_t &n, const col_t &minNodes);
void IK_GPX(const col_t &n, const col_t &minNodes, bool *R, bool *P, bool *X, col_t sizeR, col_t sizeP, col_t sizeX);
