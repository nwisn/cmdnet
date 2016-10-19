# cmdnet
 Version: 0.9 Author: Nicholas Wisniewski

The `cmdnet` package computes correlation matrix distance statistics for an ensemble of networks.

To install, first install the `devtools` package, and then type: `install_github("nwisn/cmdnet")`.

Interest in integrating multiple public datasets is growing as more scientific data is shared. Network analysis tools are important in analysis of high-dimensional data such as genomics, but difficulties arise when faced with multiple sources of data due to large variability across datasets. This algorithm, `cmdnet`, attempts to understand this variability by computing correlation matrix distance (CMD) statistics for an ensemble of networks. It can be used to quantify overall and node-specific differences between network datasets, and rank nodes by consistency of edge weights. In this way, it can also be used as a noise reduction tool by identifying highly variable nodes for removal.

