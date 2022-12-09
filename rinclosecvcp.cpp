#include "rinclosecvcp.h"

float runRInCloseCVCP(const dataset_t &D, const row_t &n, const col_t &m, const col_t &minCol)
{
	// Allocating memory for the global variable
	g_RWp = new pair<data_t, row_t>[n];


	clock_t clocks = clock();

	// Creating the supremum
	pbic_t bic = new bic_t;
	bic->A = new row_t[n];
	for (row_t i = 0; i < n; ++i)
		bic->A[i] = i;
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
	bic->biggerSup = g_biggerMinsup;

	RInCloseCVCP(D, m, minCol, bic); // call RIn-Close

	clocks = clock() - clocks;
	return ((float)clocks) / CLOCKS_PER_SEC;
}

void RInCloseCVCP(const dataset_t &D, const col_t &m, const col_t &minCol, const pbic_t &bic)
{
	stack<pbic_t> children;
	row_t *support =  new row_t[g_maxLabel];
	unsigned short label;
	row_t biggerSup;

	// Iterating across the attributes
	for (col_t j = bic->col; j < m; ++j)
	{
		if (m - j + bic->sizeB < minCol)
			break;

		if (!bic->B[j] && !bic->PN[j])
		{
			row_t qnmv = 0; // number of non-missing values
			bool allequal = true;
			for (row_t i = 0; i < bic->sizeA; ++i)
			{
				if (D[bic->A[i]][j] != MVS)
				{
					g_RWp[qnmv].first = D[bic->A[i]][j];
					g_RWp[qnmv].second = bic->A[i];
					++qnmv;
					if (D[bic->A[i]][j] != D[bic->A[0]][j])
						allequal = false;
				}
			}

			if (qnmv == bic->sizeA && allequal) // if all elements are not missing values and are the same
			{
				bic->B[j] = true; //then, add the attribute j to B[r] (incremental closure)
				++bic->sizeB;
			}
			else if (bic->biggerSup != bic->sizeA && bic->sizeA > g_smallerMinsup && qnmv >= g_smallerMinsup)
			{//bic->biggerSup==bic->sizeA means conf=100%. So, I am cutting the branch when I find a bic with conf=100%.
				bool pskipJ = true, naux1 = true, naux2 = false; // can descendants skip column j ?
				col_t fcol;

				sort(g_RWp, g_RWp + qnmv);
				
				row_t p1 = 0, p2 = 0;
				while (p1 <= qnmv - g_smallerMinsup)
				{
					biggerSup = 0;
					for (unsigned short i = 0; i < g_maxLabel; ++i) support[i] = 0; // initialize vector
					label = g_classes[g_RWp[p1].second];
					support[label] = 1;
					if (support[label] >= g_minsups[label])
					{
						if (support[label] > biggerSup) biggerSup = support[label];
					}
					while (p2 < qnmv - 1 && g_RWp[p2 + 1].first == g_RWp[p1].first)
					{
						++p2;
						label = g_classes[g_RWp[p2].second];
						support[label]++;
						if (support[label] >= g_minsups[label])
						{
							if (support[label] > biggerSup) biggerSup = support[label];
						}
					}
					if (biggerSup > 0)
					{
						pskipJ = false;
						bool iscan = RCVCP_IsCanonical(D, j, p1, p2, bic, fcol);
						if (iscan)
						{
							naux1 = false;
							pbic_t child = new bic_t;
							child->sizeA = p2 - p1 + 1;
							child->A = new row_t[child->sizeA];
							for (row_t i2 = p1; i2 <= p2; ++i2)
								child->A[i2-p1] = g_RWp[i2].second;
							child->biggerSup = biggerSup;
							child->col = j + 1;
							children.push(child);
						}
						else if (fcol >= bic->col) naux1 = false;
						else naux2 = true;
					}
					p1 = p2 + 1;
					p2 = p1;
				}
				bic->PN[j] = pskipJ || (naux1 && naux2);
			}
			else if (qnmv < g_smallerMinsup) bic->PN[j] = true; // descendants can skip column j
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
		pbic_t child = children.top();
		child->B = new bool[m];
		child->PN = new bool[m];
		for (col_t j = 0; j < m; ++j)
		{
			child->B[j] = bic->B[j];
			child->PN[j] = bic->PN[j];
		}
		child->B[child->col - 1] = true;
		child->sizeB = bic->sizeB + 1;
		RInCloseCVCP(D, m, minCol, child);
		children.pop();
	}
	delete[] bic->B;
	delete[] bic->PN;
	delete bic;

	delete[] support;
}

bool RCVCP_IsCanonical(const dataset_t &D, const col_t &y, const row_t &p1, const row_t &p2, const pbic_t &bic, col_t &fcol)
{
	row_t i;
	for (col_t j = 0; j < y; ++j)
	{
		fcol = j; // se IsCanonical=F, eh usado para guardar a coluna de falha
		if (!bic->B[j] && !bic->PN[j])
		{
			if (p1 != p2)
			{
				if (D[g_RWp[p1].second][j] != MVS)
				{
					for (i = p1 + 1; i <= p2; ++i)
					{
						if (D[g_RWp[i].second][j] == MVS || D[g_RWp[i].second][j] != D[g_RWp[p1].second][j])
							break;
					}
					if (i > p2)
						return false;
				}
			}
			else if (D[g_RWp[p1].second][j] != MVS)
				return false;
		}
	}
	return true;
}
