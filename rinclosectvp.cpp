#include "rinclosectvp.h"

float runRInCloseCTVP(const dataset_t &D, const row_t &n, const col_t &m, const row_t &minRow, const col_t &minCol)
{
	// Allocating memory for the global variable
	g_RW[0] = new row_t[n];


	clock_t clocks = clock();

	// Get the distinct values in the dataset
	set<data_t> dvalues;
	for (row_t i = 0; i < n; ++i)
	{
		for (col_t j = 0; j < m; ++j)
		{
			if (D[i][j] != MVS) dvalues.insert(D[i][j]);
		}
	}
	
	// Creates a supremum for each distinct value in the dataset
	set<data_t>::iterator it;	
	for (it=dvalues.begin(); it!=dvalues.end(); ++it)
	{
		pbic_t bic = new bic_t;
		bic->A = new row_t[n];
		for (row_t i = 0; i < n; ++i) bic->A[i] = i;
		bic->sizeA = n;
		bic->B = new bool[m];
		bic->PN = new bool[m];
		for (col_t i = 0; i < m; ++i)
		{
			bic->B[i] = false;
			bic->PN[i] = false;
		}
		bic->sizeB = 0;
		bic->col = 0;
		bic->value = *it;

		RInCloseCTVP(D, m, minRow, minCol, bic); // call RIn-Close	
	}

	clocks = clock() - clocks;
	return ((float)clocks) / CLOCKS_PER_SEC;
}

void RInCloseCTVP(const dataset_t &D, const col_t &m, const row_t &minRow, const col_t &minCol, const pbic_t &bic)
{
	queue<pbic_t> children;

	// Iterating across the attributes
	for (col_t j = bic->col; j < m; ++j)
	{
		if (m - j + bic->sizeB < minCol)
			break;

		if (!bic->B[j] && !bic->PN[j])
		{
			// Computing RW
			row_t sizeRW = 0;
			for (row_t i = 0; i < bic->sizeA; ++i)
			{
				if (D[bic->A[i]][j] == bic->value)
					g_RW[0][sizeRW++] = bic->A[i];
			}

			if (sizeRW == bic->sizeA)
			{
				bic->B[j] = true;
				++bic->sizeB;
			}
			else
			{
				if (sizeRW >= minRow)
				{
					col_t fcol;
					if (RCTVP_IsCanonical(D, j, sizeRW, bic, fcol))
					{
						pbic_t child = new bic_t;
						child->sizeA = sizeRW;
						child->A = new row_t[sizeRW];
						for (row_t i = 0; i < sizeRW; ++i)
							child->A[i] = g_RW[0][i];
						child->col = j + 1;
						child->value = bic->value;
						children.push(child);
					}
					else if (fcol < bic->col) bic->PN[j] = true;
				}
				else bic->PN[j] = true;
			}			
		}
	}

	// For the experiments about memory RAM usage
	struct rusage r_usage;
    getrusage(RUSAGE_SELF,&r_usage);
	if (r_usage.ru_maxrss > g_ram) g_ram = r_usage.ru_maxrss;

	// imprime o bicluster e elimina suas linhas da memoria
	if (bic->sizeB >= minCol)
		printBic(bic, m);
	delete[] bic->A;

	// fechando os filhos
	while (!children.empty())
	{
		pbic_t child = children.front();
		child->B = new bool[m];
		child->PN = new bool[m];
		for (col_t j = 0; j < m; ++j)
		{
			child->B[j] = bic->B[j];
			child->PN[j] = bic->PN[j];
		}
		child->B[child->col - 1] = true;
		child->sizeB = bic->sizeB + 1;
		RInCloseCTVP(D, m, minRow, minCol, child);
		children.pop();
	}
	delete[] bic->B;
	delete[] bic->PN;
	delete bic;
}

bool RCTVP_IsCanonical(const dataset_t &D, const col_t &y, const row_t &sizeRW, const pbic_t &bic, col_t &fcol)
{
	row_t i;
	for (col_t j = 0; j < y; ++j)
	{
		if (!bic->B[j] && !bic->PN[j])
		{
			for (i = 0; i < sizeRW; ++i)
				if (D[g_RW[0][i]][j] != bic->value)
					break;
			fcol = j; // if IsCanonical==false, it is used to keep the column where the test fail
			if (i == sizeRW)
				return false;
		}
	}
	return true;
}
