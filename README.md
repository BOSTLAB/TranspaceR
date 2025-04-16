# TranspaceR

TranspaceR is an R-package dedicated to process various spatial data to establish an atlas of Cellular Spatial Structures. This project focuses on two primary analyses:

1. **KNN-Based Clustering**: We employ a K-Nearest Neighbors (KNN) approach to cluster cells based on their spatial transcriptomic profiles.
2. **Cell Annotation**: Cell types are annotated using the Scimilarity tool, which documentation's can be found [here](https://genentech.github.io/scimilarity/index.html).

## Installation

TranspaceR can be installed from the source file :

```R
devtools::install_local("Path/to/TranspaceR_1.0.0.tar.gz",dependencies = T)
```

It can also be installed using devtools :
`R
devtools::install_github("BOSTLAB/TranspaceR")
`
