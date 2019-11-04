# RIn-Close

Enumerative algorithms for biclustering.

Related works:
- Veroneze, R., & Von Zuben, F. J. (2018). RIn-Close_CVC2: an even more efficient enumerative algorithm for biclustering of numerical datasets. arXiv preprint arXiv:1810.07725.
- Veroneze, R., & Von Zuben, F. J. (2017). Efficient mining of maximal biclusters in mixed-attribute datasets. arXiv preprint arXiv:1710.03289.
- Veroneze, R., Banerjee, A., & Von Zuben, F. J. (2017). Enumerating all maximal biclusters in numerical datasets. Information Sciences, 379, 288-309.
- Veroneze, R. (2016). Enumerating all maximal biclusters in numerical datasets= Enumerando todos os biclusters maximais em conjuntos de dados numéricos. PhD Thesis.



## Compiling
To compile the program, run the file  './MakeFile.sh' in the directory containing the source files.


## Running
To run the program, type './RInClose' and the arguments:

1 - Dataset's filename;

2 - Algorithm (options: cvcp, cvc, cvcma, chvp, chvpm, chv, opsm);

3 - Minimum number of rows (minRow);

4 - Minimum number of column (minCol);

5 - Epsilon or filename for file with epsilons;

6 - Output filename for the list of biclusters;

7 - Class labels' filename (optional);

8 - Confidence [0,1] (when using class labels);

**Examples are availabe in the folder 'examples'. The example dataset 'dataset' illustrates the input format.**

###Algorithm options:

**cvcp**: perfect biclusters with constant values on columns;

**cvc**: biclusters with constant (similar) values on columns - it uses the same value of epsilon for all columns;

**cvcma**: biclusters with constant (similar) values on columns - it can use different values of epsilon accross the columns;

**chvp**: perfect biclusters with coherent values - additive model.

**chvpm**: perfect biclusters with coherent values - multiplicative model.

**chv**: biclusters with coherent values - multiplicative model.

**opsm**: order-preserving submatrices