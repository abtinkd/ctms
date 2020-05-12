## CTMS (Closed Triple Micro Structure)

This project addresses the link sign prediction problem in signed networks. In these networks connections(edges) between users(nodes) are labelled with +1 and -1. These labels can signify opposite attitudes such as trust/distrust or friendship/enmity. This project can be used as a C++ library or standalone tool to predict signs of unknown links in an arbitrary network. It has been tested on Windows 10 with Visual Studio 2019 and Linux with gcc 4.8.5.

*CTMS* is simply a graph of three nodes and related edges in a way that every two nodes are connected. Edges in these small graphs can be signed or unsigned, directed or undirected. Two CTMSs are said to be equivalent if they represent the same graph when we ignore node labels. According to [\[1\]](#acknowledgments), we can have 96 distinct types of signed directed CTMSs. For example, +o|+o|+o and o+|o+|o+ are showing the same CTMS but they are different from +o|o+|+o. In this paper, it is shown that by looking at the distribution of these 96 classes of CTMSs throughout the network, we can predict the sign of a newly formed edge.

This project provides 8(+1) sign prediction methods that are used in [\[1\]](#acknowledgments). These methods are implemented as part of the `ctms` library. The datasets that are used in this project are [Epinions](http://snap.stanford.edu/data/soc-sign-epinions.html), [Slashdot](http://snap.stanford.edu/data/soc-sign-Slashdot090221.html) and [Wikipedia](http://snap.stanford.edu/data/wiki-Elec.html) and are taken from [Stanford Dataset Collection](http://snap.stanford.edu/data/index.html#signnets). In addition, there are three main sub-projects that provide the capability of applying these methods on an arbitrary *network.txt* file with the following format:
#### Example of a network with two +,- edges
```
#\tFromNodeId\tToNodeId\tSign
3\t30\t1
3\t40\t-1
```

## Compiling CTMS

In order to compile the source codes and build `convertmain`, `extractormain`, `predictormain`, `analyzermain`
and `ctmsmain`, use following commands[*]:
```
~/snap/contrib/iasbs_sharif-ctms$ make all
```
or build them separately:
```
~/snap/contrib/iasbs_sharif-ctms/converter$ make
~/snap/contrib/iasbs_sharif-ctms/extractor$ make
~/snap/contrib/iasbs_sharif-ctms/predictor$ make
~/snap/contrib/iasbs_sharif-ctms/analyzer$ make
~/snap/contrib/iasbs_sharif-ctms/ctms$ make
```
[*] Visual studio files are also included.
## Project Structure
### [ctms](ctms)
- `ctmsmain` is an all-inclusive module to apply all of the methods in cross validation based experiments.
- `ctmsnet` is the extended version of `snap-exp/signnet` augmented with various classes of CTMS and some improvements on signnet.
- `sign_prediction` has the implementations of cross validation based experiments.
- `ml` has the implementations for logistic regression optimizer.

### [extractor](extractor)
- This project can be used to generate *train* and *test* datasets from an arbitrary network. The arguments for `extractormain`, specifiy the properties of the edges that are considered for the *test* dataset. Rest of the edges in the network are used for the *train* dataset.
- Arguments:
```
arg1 -- source file path
arg2 -- output file name
arg3 -- minimum embeddedness
arg4 -- balanced or not?(1/0) [optional]
arg5 -- 0< sample size <=100 [optional]
arg6 -- maximum embeddedness [optional]
```

### [predictor](predictor)
- This project can be used to predict the signs of unlabeled edges (test data) using methods discussed in the paper. It can also be used as a guide to import the library in other projects.
- Arguments:
```
arg1 -- algorithm:
   1 <- generative-based
   2 <- receptive-based
   3 <- compound-based
   4 <- weighted generative receptive combination
   5 <- heuristic balance
   6 <- heuristic status
   7 <- logistic regression
   8 <- CTMS
   9 <- Local CTMS (not included in the paper-- SLOW!!!)
arg2 -- train file path
arg3 -- test file path
arg4 -- predictions file name
```

### [analyzer](analyzer)
- This project can be used to benchmark the performance of each method based on true labels and predicted ones.
- Arguments:
```
arg1 -- true lables file path
arg2 -- prediction file path
```

### [converter](converter)
- `convertmain` and `convert_wikipedia` are used to convert wikipedia file format.

## Acknowledgments

If you have used this software, please cite following paper:

**[1] Khodadadi, A., & Jalili, M. (2017). Sign prediction in social networks based
on tendency rate of equivalent micro-structures. Neurocomputing, 257, 175-184.**


## Contacts


For software sources, databases, helps and bugs please send me an [email](mailto:abt.kod@gmail.com)
