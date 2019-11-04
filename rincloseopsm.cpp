#include "stdafx.h"
#include "rincloseopsm.h"

float runRInCloseOPSM(const dataset_t &D, const row_t &n, const col_t &m, const row_t &minRow, const col_t &minCol)
{
	// Allocating memory for the adjacency matrix
	g_Adj = new bool*[m];
	for (col_t j = 0; j < m; ++j)
		g_Adj[j] = new bool[m];
	
	
	clock_t clocks;
	col_t minCol2 = minCol*(minCol - 1) / 2;

	// --- Step 1 - Compute the augmented matrix ---------------------
	clocks = clock();
	datasetOP_t Da;
	col_t mAug = m *(m - 1) / 2, col;
	Da = new dataOP_t*[n];
	for (row_t i = 0; i < n; ++i)
	{
		Da[i] = new dataOP_t[mAug];
		for (col_t j1 = 0; j1 < m - 1; ++j1)
		{
			for (col_t j2 = j1 + 1; j2 < m; ++j2)
			{
				col = S3_getColDa(m, j1, j2);
				
				if (D[i][j1] != MVS && D[i][j2] != MVS)
				{
					if (D[i][j1] < D[i][j2])
						Da[i][col] = 1;
					else if (D[i][j1] == D[i][j2])
						Da[i][col] = 2;
					else
						Da[i][col] = 3;
				}
				else
					Da[i][col] = MVS;
			}
		}
	}
	// -----------------------------------------------------------------

	// --- Step 2 e 3 ---------------------------------------------------
	runRInCloseCVCP_OP(Da, n, mAug, minRow, minCol2);
	// ------------------------------------------------------------------

	clocks = clock() - clocks;
	return ((float)clocks) / CLOCKS_PER_SEC;
}
