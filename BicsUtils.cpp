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
	if (g_minConf > 0 && bic->biggerSup/(double)bic->sizeA < g_minConf)
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

void closePrintFile()
{
	g_filebics.close();
}
