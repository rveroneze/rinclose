#include "stdafx.h"
#include "globalsv.h"
#include "rinclosecvcp.h"
#include "rinclosechvp.h"
#include "rinclosecvc.h"
#include "rinclosechv.h"
#include "rinclosechvpm.h"
#include "rincloseopsm.h"
#include "string.h"

bool readDataset(const string &dataSetName, dataset_t &matrix, row_t &n, col_t &m);
void printData(const dataset_t &matrix, const row_t &n, const col_t &m);
bool readEpsilons(const string &filename, data_t *epsilons, const col_t &m);
bool readAttrWeights(const string &filename, long *attrWeight, const col_t &m);
bool readClassLabels(const string &fileName, const row_t &n);
bool readConfigFile();

int main(int argc, char* argv[])
{
	if (argc < 9 || argc == 10 || argc == 11 || argc > 12)
	{
		cout << "\n!!! Wrong Arguments !!!" << endl << endl;
		cout << "List of the arguments:" << endl;
		cout << "1 - Dataset's filename;" << endl;
		cout << "2 - Algorithm (options: cvcp, cvc, cvcma, chvp, chvpm, chv, opsm);" << endl;
		cout << "3 - minRow;" << endl;
		cout << "4 - minCol;" << endl;
		cout << "5 - Epsilon or filename for file with epsilons;" << endl;
		cout << "6 - Output filename for the list of biclusters;" << endl;
		cout << "7 - minAttrWeight;" << endl;
		cout << "8 - Filename for the list of attribute weights;" << endl;
		cout << "9 - Class labels' filename (optional);" << endl;
		cout << "10 - Confidence [0,1] (when using class labels);" << endl;
		cout << "11 - Ignore biclusters with label x = ? (when using class labels);" << endl;
		exit(1);
	}

	if (!readConfigFile())
	{
		cout << "\nConfiguration file was not loaded.\n";
		cout << "Default configuration:" << endl;
		cout << "Value that represents a Missing Value: " << MVS << endl;
		cout << "Output format (1 - matlab; 2 - python): " << g_output << endl;
	}

	row_t minRow = atoi(argv[3]);
	col_t minCol = atoi(argv[4]);
	data_t epsilon, *epsilons = NULL;
	long minAttrWeight = atol(argv[7]);

	// List the user parameters
	cout << "\nArguments: " << endl;
	cout << "Dataset's filename: " << argv[1] << endl;
	cout << "Algorithm: " << argv[2] << endl;
	cout << "minRow: " << minRow << endl;
	cout << "minCol: " << minCol << endl;
	cout << "Epsilon: " << argv[5] << endl;
	cout << "File with the list of bicluster: " << argv[6] << endl;
	cout << "minAttrWeight: " << minAttrWeight << endl;
	cout << "File with the list of attribute weights: " << argv[8] << endl;
	if (argc > 9)
	{
		cout << "Class labels' filename: " << argv[9] << endl;
		cout << "Confidence: " << argv[10] << endl;
		cout << "Ignore biclusters with label x = "  << argv[11] << endl;
	}

	dataset_t matrix; // pointer to the dataset
	row_t n; // number of dataset's rows
	col_t m; // number of dataset's columns
	if (!readDataset(argv[1], matrix, n, m))
	{
		cout << "\nDataset was not loaded!";
		exit(1);
	}
	printf("\nDataset loaded: %dx%d\n\n", n, m);
	//printData(matrix, n, m);

	if (strcmp(argv[2], "cvcma") != 0)
		epsilon = atof(argv[5]);
	else
	{
		epsilons = new data_t[m];
		if (!readEpsilons(argv[5], epsilons, m))
		{
			cout << "File with the epsilons was not loaded!";
			exit(1);
		}
		cout << "File with the epsilons loaded\n\n";
	}

	long *attrWeights = new long[m];
	if (!readAttrWeights(argv[8], attrWeights, m))
	{
		cout << "File with the list of attribute weights was not loaded!";
		exit(1);
	}
	cout << "File with the list of attribute weights loaded\n\n";

	if (argc > 9)
	{
		// Le as classes dos objetos
		g_classes = new unsigned short[n];
		if (!readClassLabels(argv[9], n))
		{
			cout << "\nClass labels' file was not loaded!";
			exit(1);
		}
		printf("Class labels loaded\n\n");
		
		g_minConf = atof(argv[10]);
		g_ignoreLabel = atoi(argv[11]);
	}

	float tempo;
	openPrintFile(argv[6]);
	cout << "\nRunning..." << endl;
	if (strcmp(argv[2], "cvcp") == 0)
		tempo = runRInCloseCVCP(matrix, n, m, minRow, minCol, minAttrWeight, attrWeights);
	else if (strcmp(argv[2], "chvp") == 0)
		tempo = runRInCloseCHVP(matrix, n, m, minRow, minCol);
	else if (strcmp(argv[2], "cvc") == 0)
		tempo = runRInCloseCVC(matrix, n, m, minRow, minCol, epsilon, false);
	else if (strcmp(argv[2], "cvcma") == 0)
		tempo = runRInCloseCVCve(matrix, n, m, minRow, minCol, epsilons, minAttrWeight, attrWeights);
	else if (strcmp(argv[2], "chv") == 0)
		tempo = runRInCloseCHV(matrix, n, m, minRow, minCol, epsilon);
	else if (strcmp(argv[2], "chvpm") == 0)
		tempo = runRInCloseCHVPM(matrix, n, m, minRow, minCol);
	else if (strcmp(argv[2], "opsm") == 0)
		tempo = runRInCloseOPSM(matrix, n, m, minRow, minCol);
	closePrintFile();

	cout << "\n\nRuntime(s): " << tempo << endl;
	cout << "Number of biclusters = " << g_cont << endl;
	cout << "Memory RAM usage = " << g_ram << endl;

	// Guarda qtd de bics e runtime num txt
	// Tambem guarda ram maxima usada (em kbytes) para fins experimentais
	ofstream myfile;
	myfile.open("summary.txt", ofstream::app);
	myfile << g_cont << '\t';
	myfile << tempo << '\t';
	myfile << g_ram << endl;
	myfile.close();

	//system("pause");
	return 0;
}

bool readDataset(const string &dataSetName, dataset_t &matrix, row_t &n, col_t &m)
{
	ifstream myStream;
	myStream.open(dataSetName, ifstream::in);

	if (!myStream.is_open())
		return false;

	//Discovering the number of rows
	n = count(istreambuf_iterator<char>(myStream), istreambuf_iterator<char>(), '\n');

	//Discovering the number of columns
	data_t dbltmp;
	string line;
	m = 0;
	myStream.seekg(0);
	getline(myStream, line);
	stringstream stream(line);
	while (stream.good())
	{
		stream >> dbltmp;
		++m;
	}

	//Allocating memory
	matrix = new data_t*[n];
	for (row_t i = 0; i < n; ++i)
		matrix[i] = new data_t[m];

	//Storing the data
	myStream.seekg(0);
	for (row_t i = 0; i < n; ++i)
	{
		for (col_t j = 0; j < m; ++j)
			myStream >> matrix[i][j];
	}

	myStream.close();
	return true;
}

void printData(const dataset_t &matrix, const row_t &n, const col_t &m)
{
	for (row_t i = 0; i < n; ++i)
	{
		for (col_t j = 0; j < m; ++j)
			cout << matrix[i][j] << '\t';
		cout << endl;
	}
}

bool readEpsilons(const string &filename, data_t *epsilons, const col_t &m)
{
	ifstream myStream;
	myStream.open(filename, ifstream::in);

	if (!myStream.is_open())
		return false;

	//Storing the data
	for (col_t j = 0; j < m; ++j)
		myStream >> epsilons[j];

	myStream.close();
	return true;
}

bool readAttrWeights(const string &filename, long *attrWeights, const col_t &m)
{
	ifstream myStream;
	myStream.open(filename, ifstream::in);

	if (!myStream.is_open())
		return false;
	
	//Storing the data
	for (col_t j = 0; j < m; ++j)
		myStream >> attrWeights[j];
	
	myStream.close();
	return true;
}

bool readClassLabels(const string &fileName, const row_t &n)
{
	// Read tha class label of each object, and
	// set g_maxLabel

	g_maxLabel = 0;

	ifstream myStream;
	myStream.open(fileName, ifstream::in);

	if (!myStream.is_open())
		return false;

	//Storing the data
	myStream.seekg(0);
	for (row_t i = 0; i < n; ++i)
	{
		myStream >> g_classes[i];
		if (g_classes[i] > g_maxLabel) g_maxLabel = g_classes[i];
	}

	myStream.close();
	++g_maxLabel;

	return true;
}

bool readConfigFile()
{
	unordered_map<string,string> params;
	ifstream myStream;
	string line;

	myStream.open("config.txt", ifstream::in);

	if (!myStream.is_open())
		return false;
		

	myStream.seekg(0);
	while (myStream.good())
	{
		getline(myStream, line);
		size_t found = line.find("=");
		if (found!=string::npos)
		{
			params[line.substr(0,found)] = line.substr(found+1);
		}
	}

	cout << "myparams contains:" << endl;
	for (auto& x: params)
		cout << x.first << ": " << x.second << endl;
	
	// Setando as variaveis com os parametros
	MVS = stoi(params["MVS"]);
	if (params["OUT"].compare("matlab") == 0)
		g_output = 1;
	else
		g_output = 2;

	myStream.close();
	return true;
}