#include "stdafx.h"
#include "BicsUtils.h"

void openPrintFile(const string &filename)
{
	g_filebics.open(filename);
	if (g_output == 2) // python
	{
		g_filebics << "#!/usr/bin/env python3" << endl << endl;
		g_filebics << "bics = []" << endl;
	}
}

void printBic(const pbic_t &bic, const col_t m)
{
	if (g_minConf > 0 && getMinConf(bic->A, bic->sizeA) < g_minConf)
		return;
	
	++g_cont;
	if (g_output == 1) // matlab
	{
		g_filebics << "A{" << g_cont << "} = [";
		for (row_t i = 0; i < bic->sizeA; ++i)
			g_filebics << bic->A[i] + 1 << " ";
		g_filebics << "];\nB{" << g_cont << "} = [";
		for (col_t i = 0; i < m; ++i)
		{
			if (bic->B[i])
				g_filebics << i + 1 << " ";
		}
		g_filebics << "];\n";
	}
	else // python
	{
		g_filebics << "bics.append([[";
		for (row_t i = 0; i < bic->sizeA; ++i)
			g_filebics << bic->A[i]  << ",";
		g_filebics << "],[";
		for (col_t i = 0; i < m; ++i)
		{
			if (bic->B[i])
				g_filebics << i << ",";
		}
		g_filebics << "]])" << endl;
	}
}

void printBic(const pbic_t &bic, const col_t *B, const col_t sizeB)
{
	if (g_minConf > 0 && getMinConf(bic->A, bic->sizeA) < g_minConf)
		return;

	++g_cont;
	if (g_output == 1) // matlab
	{
		g_filebics << "A{" << g_cont << "} = [";
		for (row_t i = 0; i < bic->sizeA; ++i)
			g_filebics << bic->A[i] + 1 << " ";
		g_filebics << "];\nB{" << g_cont << "} = [";
		for (col_t i = 0; i < sizeB; ++i)
			g_filebics << B[i] + 1 << " ";
		g_filebics << "];\n";
	}
	else // python
	{
		g_filebics << "bics.append([[";
		for (row_t i = 0; i < bic->sizeA; ++i)
			g_filebics << bic->A[i] << ",";
		g_filebics << "],[";
		for (col_t i = 0; i < sizeB; ++i)
			g_filebics << B[i] << ",";
		g_filebics << "]])" << endl;
	}
}

void printBic(const pbic_t &bic, const clique_t &clique, const col_t *B)
{
	if (g_minConf > 0 && getMinConf(bic->A, bic->sizeA) < g_minConf)
		return;

	++g_cont;
	if (g_output == 1) // matlab
	{
		g_filebics << "A{" << g_cont << "} = [";
		for (row_t i = 0; i < bic->sizeA; ++i)
			g_filebics << bic->A[i] + 1 << " ";
		g_filebics << "];\nB{" << g_cont << "} = [";
		for (col_t j = 0; j < clique.size; ++j)
			g_filebics << B[clique.nodes[j]] + 1 << " ";
		g_filebics << "];\n";
	}
	else // python
	{
		g_filebics << "bics.append([[";
		for (row_t i = 0; i < bic->sizeA; ++i)
			g_filebics << bic->A[i] << ",";
		g_filebics << "],[";
		for (col_t j = 0; j < clique.size; ++j)
			g_filebics << B[clique.nodes[j]] << ",";
		g_filebics << "]])" << endl;
	}
}


void closePrintFile()
{
	g_filebics.close();
}

double getMinConf(row_t *A, row_t size)
{
// Determina a confianca do bicluster considerando a classe que ele mais representa

	unsigned int *contClass = new unsigned int[g_maxLabel];
	for (unsigned short i = 0; i < g_maxLabel; ++i) contClass[i] = 0; // inicializa contador

	for (row_t i = 0; i < size; ++i) ++contClass[ g_classes[A[i]] ]; // contando a representatividade de cada classe
	
	unsigned int maior = 0;
	for (unsigned short i = 0; i < g_maxLabel; ++i)
	{
		if (contClass[i] > maior) maior = contClass[i];
	}

	delete [] contClass;

	return maior / (double) size;
}

