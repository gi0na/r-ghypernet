  <!-- badges: start -->
[![Generic badge](https://www.r-pkg.org/badges/version-last-release/ghypernet)](https://cran.r-project.org/package=ghypernet)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2555300.svg)](https://doi.org/10.5281/zenodo.2555300)
[![Generic badge](https://cranlogs.r-pkg.org/badges/last-month/ghypernet)](https://cran.r-project.org/package=ghypernet)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL%20v3-red.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![Travis build status](https://travis-ci.com/gi0na/r-ghypernet.svg?branch=master)](https://travis-ci.com/gi0na/r-ghypernet)
  <!-- badges: end -->
 

![](man/figures/logo.svg)

# Introduction
`ghypernet` is an OpenSource R package that allows to estimate and work with Generalised Hypergeometric Ensembles of Random Graphs (gHypEG).

`ghypernet` has been developed specifically for the analysis of networks characterised by a large number of repeated edges.
It provides efficient methods to perform hypothesis testing and model selection on such data.

Explore the provided Vignettes for some examples on how to analyse networks with `ghypernet`.

# Installation
```
# Install ghypernet from CRAN
install.packages("ghypernet")

# Or the development version from GitHub:
devtools::install_github("gi0na/r-ghypernet")
```

# Dependencies
The package uses the library `BiasedUrn` to work with Wallenius' non-central hypergeometric distribution.
Although this is not required, it is recommended to install the `BiasedUrn` R package, increasing the number of 'colors', i.e., the number of dimensions of the distribution.
It can be easily done modifying the _makevar_ file.
In case the `BiasedUrn` library cannot be found, all computations will be performed using the multinomial approximation.

# References
The theoretical foundation of the generalised hypergeometric ensemble, gHypEGs, has been developed in the following works:

Casiraghi, G., Nanumyan, V., Scholtes, I., & Schweitzer, F. (2016). [Generalized Hypergeometric Ensembles: Statistical Hypothesis Testing in Complex Networks.](https://arxiv.org/abs/1607.02441) ArXiv Preprint ArXiv:1607.02441.

Casiraghi, G. (2017). [Multiplex Network Regression: How do relations drive interactions?](https://arxiv.org/abs/1702.02048). ArXiv Preprint ArXiv:1702.02048, 15.

Casiraghi, G., Nanumyan, V., Scholtes, I., & Schweitzer, F. (2017). [From Relational Data to Graphs: Inferring Significant Links Using Generalized Hypergeometric Ensembles](https://doi.org/10.1007/978-3-319-67256-4_11) (Vol. 10540, pp. 111–120). Springer Verlag.

Casiraghi, G. (2019). [The block-constrained configuration model.](https://doi.org/10.1007/s41109-019-0241-1) Applied Network Science, 4(1), 123. 

Brandenberger, L., Casiraghi, G., Nanumyan, V., & Schweitzer, F. (2019). [Quantifying triadic closure in multi-edge social networks.](https://doi.org/10.1145/3341161.3342926) Proceedings of the 2019 IEEE/ACM International Conference on Advances in Social Networks Analysis and Mining, 307–310.

Casiraghi, G., & Nanumyan, V. (2021). [Configuration models as an urn problem.](https://www.nature.com/articles/s41598-021-92519-y) Sci Rep 11, 13416.

Casiraghi, G. (2021) [The likelihood-ratio test for multi-edge network models.](https://iopscience.iop.org/article/10.1088/2632-072X/ac0493) J. Phys. Complex. 2 035012.

# Acknowledgements
The research and development behind `ghypernet` is performed at the Chair of Systems Design, ETH Zürich.

# Contributors

[Giona Casiraghi](https://github.com/gi0na) (project lead)

[Vahan Nanumyan](https://www.sg.ethz.ch/)

[Laurence Brandenberger](https://www.sg.ethz.ch/team/people/brandenberger-laurence/)

# Copyright
`ghypernet` is licensed under the [GNU Affero General Public License](https://choosealicense.com/licenses/agpl-3.0/).

(c) Copyright ETH Zürich, 2016-2021
