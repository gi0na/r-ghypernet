[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2555300.svg)](https://doi.org/10.5281/zenodo.2555300)
[![Generic badge](https://www.r-pkg.org/badges/version/ghypernet)]()
[![Generic badge](https://cranlogs.r-pkg.org/badges/last-month/ghypernet)]()
# Introduction
`ghypernet` is an OpenSource R package that allows to estimate and work with Generalised Hypergeometric Ensembles of Random Graphs (gHypEG).

`ghypernet` has been developed specifically for the analysis of large networks characterised by a large number of repeated edges.
It provides efficient methods to perform hypothesis testing and model selection on such data.

The theoretical foundation of this paper, gHypEGs, was developed in the following works:

  1. Casiraghi, G., Nanumyan, V., Scholtes, I., Schweitzer, F. (2016) [Generalized Hypergeometric Ensembles: Statistical Hypothesis Testing in Complex Networks.](http://arxiv.org/abs/1607.02441) arXiv Prepr. arXiv1607.02441
  2. Casiraghi, G., Nanumyan, V., Scholtes, I., Schweitzer, F. (2017) [From Relational Data to Graphs: Inferring Significant Links Using Generalized Hypergeometric Ensembles](https://link.springer.com/chapter/10.1007/978-3-319-67256-4_11) in Social Informatics. SocInfo 2017 111-120 (Springer Verlag, 2017). doi:10.1007/978-3-319-67256-4_11
  3. Casiraghi, G., Nanumyan, V. (2018) [Generalised hypergeometric ensembles of random graphs: the configuration model as an urn problem.](http://arxiv.org/abs/1810.06495)  arXiv Prepr. arXiv1810.06495
  4. Casiraghi, G. (2018) [Analytical Formulation of the Block-Constrained Configuration Model.](http://arxiv.org/abs/1811.05337)  arXiv Prepr. arXiv1811.05337

# Dependencies
The package uses the library `BiasedUrn` to work with Wallenius' non-central hypergeometric distribution.
Although this is not required, it is recommended to install the `BiasedUrn` R package, increasing the number of 'colors', i.e., the number of dimensions of the distribution.
It can be easily done modifying the _makevar_ file.
In case the `BiasedUrn` library cannot be found, all computations will be performed using the multinomial approximation.

# Installation
You can install this package directly from GitHub.
In R, run the following commands to install the package:
```
install.packages('devtools')
devtools::install_github("gi0na/r-ghypernet")

library(ghypernet)
```


# Acknowledgements
The research and development behind `ghypernet` is performed at the Chair of Systems Design, ETH Zürich.

# Contributors

[Giona Casiraghi](http://giona.info) (project lead)

[Vahan Nanumyan](https://www.sg.ethz.ch/)

[Laurence Brandenberger](https://www.sg.ethz.ch/team/people/brandenberger-laurence/)

# Copyright
`ghypernet` is licensed under the [GNU Affero General Public License](https://choosealicense.com/licenses/agpl-3.0/).

(c) Copyright ETH Zürich, 2016-2020
