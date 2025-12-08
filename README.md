# RIn-Close

Enumerative algorithms for biclustering.

Related works:
- Veroneze, R., & Von Zuben, F. J. (accepted). Scalability achievements for enumerative biclustering with online partitioning: case studies involving mixed-attribute datasets. Engineering Applications of Artificial Intelligence.
- Veroneze, R., & Von Zuben, F. J. (2020). New advances in enumerative biclustering algorithms with online partitioning. arXiv preprint arXiv:2003.04726.
- Veroneze, R., & Von Zuben, F. J. (2018). RIn-Close_CVC2: an even more efficient enumerative algorithm for biclustering of numerical datasets. arXiv preprint arXiv:1810.07725.
- Veroneze, R., & Von Zuben, F. J. (2017). Efficient mining of maximal biclusters in mixed-attribute datasets. arXiv preprint arXiv:1710.03289.
- Veroneze, R., Banerjee, A., & Von Zuben, F. J. (2017). Enumerating all maximal biclusters in numerical datasets. Information Sciences, 379, 288-309.
- Veroneze, R. (2016). Enumerating all maximal biclusters in numerical datasets= Enumerando todos os biclusters maximais em conjuntos de dados numéricos. PhD Thesis.
- Veroneze, R., Banerjee, A., & Von Zuben, F. J. (2014). Enumerating all maximal biclusters in numerical datasets. arXiv preprint arXiv:1403.3562.


## Branch for testing a new feature
This branch contains an adaptation for mining **class association rules**.

The user must specify **one minimum support** value *and* **one minimum confidence** value **per class label**.

* The **minimum support** is used to prune the search space during pattern exploration.
* The **minimum confidence** acts only as a filter when deciding which rules to output.

For example, consider a dataset with two class labels, `c1` and `c2`.
For each pattern `X`, we can derive the rules `X → c1` and `X → c2`.
However, only the rule achieving the **highest confidence** is evaluated. It is retained and printed **only if** its confidence is greater than or equal to the corresponding class-specific minimum confidence threshold.


## Compiling
To compile the program, run the file  './MakeFile.sh' in the directory containing the source files.

**Examples are availabe in the folder 'examples'. The example dataset 'dataset' illustrates the input format.**
