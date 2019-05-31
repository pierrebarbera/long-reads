## Data and scripts for _"Long metabarcoding of the eukaryotic rDNA operon to phylogenetically and taxonomically resolve environmental diversity"_ (Jamy et al. 2019)

### Structure

This repository holds the data and part of the analysis stack of the above mentioned paper. It is structured in the following way:

* **[data](data/)** holds the MSAs for various combination of query sequences (either full or constrained to V4) and reference sequences, as well as taxonomic information on reference taxa, and outgroup information
* **[trees](trees/)** contains the results of tree searches based on the aforementioned alignments
* **[jplace](jplace/)** contains the result of phylogenetic placement, the first sub-dir corresponding to the reference tree used, and the second corresponding to which query sequence MSA was used
* **[assign](assign/)** contains taxonomic assignment results based either on the tree inferred from mixing references with queries (`comprehensive_*`) or on phylogenetic placement (the rest)
* **[leave1out](leave1out/)** contains the result of the leave-one-out tests
* **[visualization](visualization/)** code and results for visualizations
* **[src](src/)** source files used in analysis
* **[preliminaries](preliminaries/)** data from preliminary study of the data, included soely for transparency reasons

Additionally, the scripts to perform analysis and produce the results are located in the base directory.

### Nomenclature
The `V4` in this context means usage of the query alignments masked to only include the V4 hypervariable region.

### Citation
```
Long metabarcoding of the eukaryotic rDNA operon to phylogenetically and taxonomically resolve environmental diversity
Mahwash Jamy, Rachel Foster, Pierre Barbera, Lucas Czech, Alexey Kozlov, Alexandros Stamatakis, David Bass, Fabien Burki
bioRxiv 627828; doi: https://doi.org/10.1101/627828
```