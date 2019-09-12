## CTMS (Closed Triple Micro Structure)

*CTMS* is simply a graph of three nodes and related edges in a way that every two nodes are connected. Edges in these small graphs can be signed or unsigned, directed or undirected. Two CTMSs are said to be equivalent if they represent the same graph when we ignore node labels. According to [\[1\]](#acknowledgments), we can have 96 distinct types of signed directed CTMSs. For example, +o|+o|+o and o+|o+|o+ are showing the same CTMS but they are different from +o|o+|+o. In this paper, it is shown that by looking at the distribution of these 96 classes of CTMSs throughout the network, we can predict the sign of a newly formed edge.

This project is providing 8 sign prediction methods that are used in [\[1\]](#acknowledgments) and the libraries that are used to extract features, build models and make predictions. The datasets that are used in this project are [Epinions](http://snap.stanford.edu/data/soc-sign-epinions.html), [Slashdot](http://snap.stanford.edu/data/soc-sign-Slashdot090221.html) and [Wikipedia](http://snap.stanford.edu/data/wiki-Elec.html) and are taken from [snap-datasets](http://snap.stanford.edu/data/index.html#signnets). However, any *network.txt* file with the following format can be used as the network.
#### A Network Example with 2 +,- Edges
```
#\tFromNodeId\tToNodeId\tSign
3\t30\t1
3\t40\t-1
```

## Compiling CTMS

In order to compile the source codes and build executables for `convert_wikipedia`
and `ctmsmain`, use following command:
```
~/snap/contrib/iasbs_sharif-ctms$ make all
```
or make them separately:
```
~/snap/contrib/iasbs_sharif-ctms/convert_wikipedia$ make
~/snap/contrib/iasbs_sharif-ctms/ctms$ make
```

## Project Structure
### [ctms](ctms)
- ctmsmain:
- ctmsnet:
- sign_prediction:
- ml:

### [convert_wikipedia](convert_wikipedia)

- convertmain:
- convert_wikipedia:


## Execute CTMS example
To run project on *testnet* dataset
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

If you have used CTMS software, please cite the following paper:

**[1] Khodadadi, A., & Jalili, M. (2017). Sign prediction in social networks based
on tendency rate of equivalent micro-structures. Neurocomputing, 257, 175-184.**


## Contacts


For software sources, databases, helps and bugs please send me an [email](khodadaa[at]oregonstate[dot]edu)
