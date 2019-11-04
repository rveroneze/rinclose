#include "stdafx.h"
#include "IK_GPX.h"

void runIK_GPX(const col_t &n, const col_t &minNodes)
{
	bool *R = new bool[n], *P = new bool[n], *X = new bool[n];
	col_t sizeR = 0, sizeP = n, sizeX = 0;

	for (col_t i = 0; i < n; ++i)
	{
		R[i] = false;
		P[i] = true;
		X[i] = false;
	}

	IK_GPX(n, minNodes, R, P, X, sizeR, sizeP, sizeX);

	delete[] R;
	delete[] P;
	delete[] X;
}

void IK_GPX(const col_t &n, const col_t &minNodes, bool *R, bool *P, bool *X, col_t sizeR, col_t sizeP, col_t sizeX)
{
	if (sizeP == 0 && sizeX == 0 && sizeR >= minNodes)
	{
		clique_t no;
		no.nodes = new col_t[sizeR];
		col_t p = 0;
		for (col_t i = 0; i < n; ++i)
			if (R[i])
				no.nodes[p++] = i;
		no.size = sizeR;
		g_cliques.push(no);
	}
	else
	{
		// Finding up: the pivot vertex
		col_t pmaior = n + 1; // inicializando com qq coisa invalida
		int maior = -1, card;
		for (col_t up = 0; up < n; ++up)
		{
			if (P[up] || X[up])
			{
				card = 0;
				for (col_t j = 0; j < n; ++j)
					if (P[j] && g_Adj[up][j]) ++card;
				if (card > maior)
				{
					maior = card;
					pmaior = up;
				}
			}
		}

		for (col_t ui = 0; ui < n; ++ui)
		{
			if (P[ui] && !g_Adj[ui][pmaior])
			{
				bool *Rnew = new bool[n], *Pnew = new bool[n], *Xnew = new bool[n];
				col_t sizePnew = 0, sizeXnew = 0;

				P[ui] = false;
				--sizeP;

				for (col_t j = 0; j < n; ++j)
				{
					Pnew[j] = false;
					Xnew[j] = false;
					Rnew[j] = R[j];
					if (P[j] && g_Adj[ui][j])
					{
						Pnew[j] = true;
						++sizePnew;
					}
					if (X[j] && g_Adj[ui][j])
					{
						Xnew[j] = true;
						++sizeXnew;
					}

				}
				Rnew[ui] = true;

				if (sizeR + 1 + sizePnew >= minNodes)
					IK_GPX(n, minNodes, Rnew, Pnew, Xnew, sizeR + 1, sizePnew, sizeXnew);

				delete[] Rnew;
				delete[] Pnew;
				delete[] Xnew;

				X[ui] = true;
				
				if (sizeR + sizeP < minNodes)
					break;
			}
		}
	}
}