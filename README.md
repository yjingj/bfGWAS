Bayesian Functional Genome-wide Association Study (bfGWAS) 
GWAS tool based on a Bayeisan variable selection model that accounts for functional genomic information and LD.

C++ libraries used in this tool: zlib, gsl, eigen3, lapack, atlas, blas, libStatGen. For success compilation of the tool, please revise the pathes to these libraries accordingly.

Please use the version of bfGWAS/libStateGen/MemoryAllocators.h and bfGWAS/libStateGen/MemoryAllocators.cpp to compile libStatGen.a. The ones from original "https://github.com/statgen/libStatGen.git" will cause error.



