## CTMS (Closed Triple Micro Structure)

CTMS is simply a graph of three nodes and related edges in a way that every two 
nodes are connected. Edges in these small graphs can be signed or unsigned, 
directed or undirected. Two CTMSs are said to be equivalent if they represent
the same graph when we ignore node labels. According to [1], we can have 96 
distinct types of signed directed CTMSs. For example, +o|+o|+o and o+|o+|o+ are
showing the same CTMS but they are different from +o|o+|+o. In this paper, it is
shown that by looking at the distribution of these 96 classes of CTMSs throughout
the network, we can predict the sign of a newly formed edge.

This project is providing 8 sign prediction methods that are used in [1] and the
libraries that are used to extract features, build models and make predictions.
The datasets that are used in this project are [Epinions](http://snap.stanford.edu/data/soc-sign-epinions.html), [Slashdot](http://snap.stanford.edu/data/soc-sign-Slashdot090221.html) and [Wikipedia](http://snap.stanford.edu/data/wiki-Elec.html)
and are taken from [snap-datasets](http://snap.stanford.edu/data/index.html#signnets).
However, any *network.txt* file with the following format can be used as the network.
#### Network File Format
```
#\tFromNodeId\tToNodeId\tSign
```
Example:
```
3\t30\t1
3\t40\t-1
```

### [convert_wikipedia](/convert_wikipedida)
ctmsmain.cpp


## Compiling CTMS

In order to compile the source codes and build executables for convert_wikipedia
and ctmsmain, use following command:
```
snap/contrib/iasbs_sharif-ctms$ make all
```
or make them separately:
```
~/snap/contrib/iasbs_sharif-ctms/ctms$ make
~/snap/contrib/iasbs_sharif-ctms/convert_wikipedia$ make
```



## Execute CTMS example





## Acknowledgments

If you have used CTMS software, please cite the following paper:

[1] Khodadadi, A., & Jalili, M. (2017). Sign prediction in social networks based
on tendency rate of equivalent micro-structures. Neurocomputing, 257, 175-184.


## Contacts


For software sources, databases, helps and bugs please send me an [email](khodadaa@oregonstate.edu)
