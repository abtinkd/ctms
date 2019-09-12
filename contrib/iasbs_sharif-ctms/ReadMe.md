## CTMS (Closed Triple Micro Structure)

This project addresses the link sign prediction problem is signed networks which are a type of networks that connections(edges) between users(nodes) are labelled with +1 and -1 signs. These labels can signify opposite attitudes such as trust/distrust or friendship/enmity. This project can be used as a C++ library or standalone tool to predict signs of unknown links in an arbitrary network. It has been tested on Windows 10 with Visual Studio 2019 and Linux with gcc 4.8.5.

*CTMS* is simply a graph of three nodes and related edges in a way that every two nodes are connected. Edges in these small graphs can be signed or unsigned, directed or undirected. Two CTMSs are said to be equivalent if they represent the same graph when we ignore node labels. According to [\[1\]](#acknowledgments), we can have 96 distinct types of signed directed CTMSs. For example, +o|+o|+o and o+|o+|o+ are showing the same CTMS but they are different from +o|o+|+o. In this paper, it is shown that by looking at the distribution of these 96 classes of CTMSs throughout the network, we can predict the sign of a newly formed edge.

This project is providing 8 sign prediction methods that are used in [\[1\]](#acknowledgments) and the libraries that are used to extract features, build models and make predictions. The datasets that are used in this project are [Epinions](http://snap.stanford.edu/data/soc-sign-epinions.html), [Slashdot](http://snap.stanford.edu/data/soc-sign-Slashdot090221.html) and [Wikipedia](http://snap.stanford.edu/data/wiki-Elec.html) and are taken from [Stanford Dataset Collection](http://snap.stanford.edu/data/index.html#signnets). However, any *network.txt* file with the following format can be used as the network.
#### Example of a network with two +,- edges
```
#\tFromNodeId\tToNodeId\tSign
3\t30\t1
3\t40\t-1
```

## Compiling CTMS

In order to compile the source codes and build `convertmain`
and `ctmsmain`, use following commands:
```
~/snap/contrib/iasbs_sharif-ctms$ make all
```
or build them separately:
```
~/snap/contrib/iasbs_sharif-ctms/converter$ make
~/snap/contrib/iasbs_sharif-ctms/ctms$ make
```

## Project Structure
### [ctms](ctms)
- `ctmsmain` is the main file to run all of the methods.
- `ctmsnet` is the extended version of `snap-exp/signnet` augmented with various classes of CTMS and some improvements on signnet.
- `sign_prediction` has the implementations for sign prediction methods.
- `ml` has the implementations for logistic regression optimizer.

### [converter](converter)

- `convertmain` and `convert_wikipedia` are used to convert wikipedia file format.


## Execute CTMS example
Run project on the *testnet* dataset with following instructions:
- 1st argument after `ctmsmain` is the path to the network file.
- 2nd argument after `ctmsmain` shows the minimum embeddedness required for the edges that are involved in the predictions.
- 3rd argument after `ctmsmain` is used to force a subnet of equal positive and negative edges.

After network is built, you will be able to set the prediction method and setting a network name to log features.

```
~/snap/contrib/iasbs_sharif-ctms/ctms$ ctmsmain ./data/testnet.txt 1 0
input file: ./data/testnet.txt
minimum embeddedness: 1
limit to a random subnet of balanced edge signs: no
loading full network ...
extracting sub network ...
  fullnet edge count: 9011
  sub-net edge count: 7726
  positive edge ratio:  0.83
Methods: All(0) NAIVE(1) CTMS(2) LogReg(3) BALANCE(4) STATUS(5): 2
Set a net name to enable logging in <./results>. (press Enter to disable it): testnet

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Sign prediction: CTMS-based method (probabilistic inference, 96 CTMS classes)

Extracting edge-feature matrix from the network ...
100%  [0.39s]

Fold 1   --- --- --- --- --- --- --- --- --- --- --- ---  Time: 0.00s
Recreating train & test data... COMPLETED
Calculating Theta ...
99%[0.01s]

TPR ( 0.94 ), TNR ( 0.53 ).
Acc: 0.8769
...
```




## Acknowledgments

If you have used this software, please cite following paper:

**[1] Khodadadi, A., & Jalili, M. (2017). Sign prediction in social networks based
on tendency rate of equivalent micro-structures. Neurocomputing, 257, 175-184.**


## Contacts


For software sources, databases, helps and bugs please send me an [email](mailto:abt.kod@gmail.com)
