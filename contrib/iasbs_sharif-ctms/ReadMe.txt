===============================================================================
    CTMS (Closed Triple Micro Structure)
===============================================================================
CTMS is simply a graph of three nodes and related edges in way that every two 
nodes are connected. Edges in these small graphs can be signed or unsigned, 
directed or undirected. Two CTMSs are said to be equivalent if they represent
the same graph when we ignore node labels. According to [1], we can have 96 
distinct types of signed directed CTMSs. For example, +o|+o|+o and o+|o+|o+ are
showing the same CTMS but they are different from +o|o+|+o.In this paper, it is
shown that by looking at the distribution of these 96 classes of CTMSs throughout
the network, we can predict the sign of a newly formed edge.

This project is providing 8 sign prediction methods that are used in [1] and the
libraries that are used to extract features, build models and make predictions.
The datasets that are used in this project are Epinions, Slashdot and Wikipedia
and are taken from [snap-datasets](http://snap.stanford.edu/data/index.html#signnets).

http://snap.stanford.edu/data/wiki-Elec.html

ctmsmain.cpp

===============================================================================
    Compiling CTMS
===============================================================================



===============================================================================
    Execute CTMS example
===============================================================================



===============================================================================
    Acknowledgments
===============================================================================
If you have used CTMS software, please cite the following paper:

[1] Khodadadi, A., & Jalili, M. (2017). Sign prediction in social networks based
on tendency rate of equivalent micro-structures. Neurocomputing, 257, 175-184.

===============================================================================
    Contacts
===============================================================================

For software sources, databases, helps and bugs please send an email to
Abtin Khodadadi khodadaa@oregonstate.edu
