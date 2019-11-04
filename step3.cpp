#include "stdafx.h"
#include "step3.h"

void S3_step3(const dataset_t &Da, const row_t &n, const col_t &mAug, const col_t &minCol2, const data_t &epsilon, const pbic_t &bic)
{
	col_t m = (1 + sqrt(1 + 8 * mAug)) / 2, minCol = (1 + sqrt(1 + 8 * minCol2)) / 2;

	col_t sizeB2, *B2 = new col_t[m], *B2p = new col_t[m], col1, col2, j, j1, j2;
	bool *B2b = new bool[m];
	if (bic->sizeB >= minCol2)
	{
		// Compute B2
		sizeB2 = 0;
		for (j = 0; j < m; ++j)
		{
			B2b[j] = false;
			B2p[j] = 0;
		}
		for (j = 0; j < mAug; ++j)
		{
			if (bic->B[j])
			{
				S3_getColsOrig(m, j, col1, col2);
				B2b[col1] = true;
				B2b[col2] = true;
			}
		}
		for (j = 0; j < m; ++j)
		{
			if (B2b[j])
			{
				B2p[j] = sizeB2;
				B2[sizeB2] = j;
				++sizeB2;
			}
		}

		// Compute the adjacency matrix
		for (j1 = 0; j1 < sizeB2; ++j1)
		for (j2 = 0; j2 < sizeB2; ++j2)
			g_Adj[j1][j2] = false;
		for (j = 0; j < mAug; ++j)
		{
			if (bic->B[j])
			{
				S3_getColsOrig(m, j, col1, col2);
				g_Adj[B2p[col1]][B2p[col2]] = true;
				g_Adj[B2p[col2]][B2p[col1]] = true;
			}
		}

		// Find all maximal cliques
		runIK_GPX(sizeB2, minCol);

		// Compute the CHV biclusters
		while (!g_cliques.empty())
		{
			clique_t clique = g_cliques.front();
			if (sizeB2 == clique.size)
				printBic(bic, B2, sizeB2);
			else
			{
				// verifico maximalidade
				if (S3_IsMaximal(Da, n, m, epsilon, bic, clique, B2))
					printBic(bic, clique, B2);
			}
			delete[] clique.nodes;
			g_cliques.pop();
		}
	}
	delete[] B2;
	delete[] B2b;
	delete[] B2p;
}

bool S3_IsMaximal(const dataset_t &Da, const row_t &n, const col_t &m, const data_t &epsilon, const pbic_t &bic, const clique_t &clique, const col_t *B2)
{
	data_t maior, menor;
	col_t j1, j2, colda;
	row_t i, k;

	bool *rows = new bool[n];
	for (i = 0; i < n; ++i)
		rows[i] = false;
	for (i = 0; i < bic->sizeA; ++i)
		rows[bic->A[i]] = true;

	for (i = 0; i < n; ++i)
	{
		if (!rows[i])
		{
			for (j1 = 0; j1 < clique.size - 1; ++j1)
			{
				for (j2 = j1 + 1; j2 < clique.size; ++j2)
				{
					colda = S3_getColDa(m, B2[clique.nodes[j1]], B2[clique.nodes[j2]]);
					if (Da[i][colda] == MVS)
						break;
					maior = Da[i][colda];
					menor = maior;
					for (k = 0; k < bic->sizeA; ++k)
					{
						if (Da[bic->A[k]][colda] == MVS)
							break;
						if (Da[bic->A[k]][colda] > maior)
							maior = Da[bic->A[k]][colda];
						if (Da[bic->A[k]][colda] < menor)
							menor = Da[bic->A[k]][colda];
						if (maior - menor > epsilon)
							break;
					}
					if (k < bic->sizeA)
						break;
				}
				if (j2 < clique.size)
					break;
			}
			if (j1 == clique.size - 1)
				return false;
		}
	}
	return true;
}

col_t S3_getColDa(const col_t &m, const col_t &col1, const col_t &col2)
{
	return (col_t)(col1 * (m - (col1 + 1) / 2.0) + col2 - col1 - 1);
}

void S3_getColsOrig(const col_t &m, const col_t &colAug, col_t &c1, col_t &c2)
{
	c1 = (col_t)(floor(m + 0.5 - sqrt(pow(m, 2) - m + 0.25 - 2 * (colAug))) - 1);
	c2 = (col_t)(colAug - c1 * (m - (c1 + 1) / 2.0) + c1 + 1);
}

void S3_step3_OP(const datasetOP_t &Da, const row_t &n, const col_t &mAug, const col_t &minCol2, const pbic_t &bic)
{
	col_t m = (1 + sqrt(1 + 8 * mAug)) / 2, minCol = (1 + sqrt(1 + 8 * minCol2)) / 2;

	col_t sizeB2, *B2 = new col_t[m], *B2p = new col_t[m], col1, col2, j, j1, j2;
	bool *B2b = new bool[m];
	if (bic->sizeB >= minCol2)
	{
		// Compute B2
		sizeB2 = 0;
		for (j = 0; j < m; ++j)
		{
			B2b[j] = false;
			B2p[j] = 0;
		}
		for (j = 0; j < mAug; ++j)
		{
			if (bic->B[j])
			{
				S3_getColsOrig(m, j, col1, col2);
				B2b[col1] = true;
				B2b[col2] = true;
			}
		}
		for (j = 0; j < m; ++j)
		{
			if (B2b[j])
			{
				B2p[j] = sizeB2;
				B2[sizeB2] = j;
				++sizeB2;
			}
		}

		// Compute the adjacency matrix
		for (j1 = 0; j1 < sizeB2; ++j1)
		for (j2 = 0; j2 < sizeB2; ++j2)
			g_Adj[j1][j2] = false;
		for (j = 0; j < mAug; ++j)
		{
			if (bic->B[j])
			{
				S3_getColsOrig(m, j, col1, col2);
				g_Adj[B2p[col1]][B2p[col2]] = true;
				g_Adj[B2p[col2]][B2p[col1]] = true;
			}
		}

		// Find all maximal cliques
		runIK_GPX(sizeB2, minCol);

		// Compute the OPSM biclusters
		while (!g_cliques.empty())
		{
			clique_t clique = g_cliques.front();
			if (sizeB2 == clique.size)
				printBic(bic, B2, sizeB2);
			else
			{
				// verifico maximalidade
				if (S3_IsMaximal_OP(Da, n, m, bic, clique, B2))
					printBic(bic, clique, B2);
			}
			delete[] clique.nodes;
			g_cliques.pop();
		}
	}
	delete[] B2;
	delete[] B2b;
	delete[] B2p;
}

bool S3_IsMaximal_OP(const datasetOP_t &Da, const row_t &n, const col_t &m, const pbic_t &bic, const clique_t &clique, const col_t *B2)
{
	data_t maior, menor;
	col_t j1, j2, colda;
	row_t i, k;

	bool *rows = new bool[n];
	for (i = 0; i < n; ++i)
		rows[i] = false;
	for (i = 0; i < bic->sizeA; ++i)
		rows[bic->A[i]] = true;

	for (i = 0; i < n; ++i)
	{
		if (!rows[i])
		{
			for (j1 = 0; j1 < clique.size - 1; ++j1)
			{
				for (j2 = j1 + 1; j2 < clique.size; ++j2)
				{
					colda = S3_getColDa(m, B2[clique.nodes[j1]], B2[clique.nodes[j2]]);
					if (Da[i][colda] == MVS)
						break;
					maior = Da[i][colda];
					menor = maior;
					for (k = 0; k < bic->sizeA; ++k)
					{
						if (Da[bic->A[k]][colda] == MVS)
							break;
						if (Da[bic->A[k]][colda] > maior)
							maior = Da[bic->A[k]][colda];
						if (Da[bic->A[k]][colda] < menor)
							menor = Da[bic->A[k]][colda];
						if (maior - menor > 1)
							break;
					}
					if (k < bic->sizeA)
						break;
				}
				if (j2 < clique.size)
					break;
			}
			if (j1 == clique.size - 1)
				return false;
		}
	}
	return true;
}
