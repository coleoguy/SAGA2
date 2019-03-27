SAGA2
====

Software for the Analysis of Genetic Architecture

The package SAGA was a collection functions to ease and hopefully improve the quality of line cross analysis of genetic architecture.  The overall goal was to allow for an easy and straightforward implementation of model averaged analysis composite genetic effects using line means.  SAGA2 takes the next step in the analysis of this type of data by building in environemntal variation and making the actual analysis easier for the end user.

In the past users had to make or confirm the type of C-matrix that correctly described the composite genetic effects in their experiment.  With SAGA2 users simply supply the breeding design of their experiment and the software generates the appropriate C-matrix automatically.

These insrtuction will install SAGA2: 

`install.packages("devtools")`

`library(devtools)`

`install_github('coleoguy/SAGA2')`

`library(SAGA2)`

If you have questions or problems please let me know [coleoguy@gmail.com](mailto:coleoguy@gmail.com).
