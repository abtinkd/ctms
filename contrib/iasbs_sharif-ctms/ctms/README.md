## Example for cross validation based experiments
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
