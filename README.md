
## About

As more scientific data is shared, better tools are needed for integrating multiple public datasets. Network analysis tools are important in analysis of high-dimensional data such as genomics, but difficulties arise when faced with multiple sources of data due to large variability across datasets. This algorithm, `cmdnet`, attempts to understand this variability by computing correlation matrix distance (CMD) statistics for an ensemble of networks. It can be used to quantify overall and node-specific differences between network datasets, and rank nodes by consistency of edge weights. In this way, it can also be used as a noise reduction tool by identifying highly variable nodes for removal.

## Vignette
There is a [vignette](cmdnet_vignette.pdf) illustrating the use of `cmdnet`.

## Installation
To install, first install the `devtools` package, and then type:
`install_github("nwisn/cmdnet")`.



