# DecayFitting
Code to simulate beta (and beta-n + beta-2n) decay of isomeric chain & fit those components

The simulation code is MCGenRoot.C and should be run within a ROOT environment. There are many parameters that require tuning for your specific case, and it should be easy enough to follow to add more. This assumes we start with nucleus A which can decay via beta, beta-n, and beta-2n decays. The naming convention is:

A - parent

B - 1beta-0n daughter

C - 1beta-1n daughter

D - 1beta-2n daughter

E - 2beta-0n granddaughter

F - 2beta-1n granddaughter

G - 2beta-2n granddaughter

H - 2beta-3n granddaughter

I - 3beta-0n great-granddaughter

J - 3beta-1n great-granddaughter

K - 3beta-2n great-granddaughter

L - 3beta-3n great-granddaughter

M - 4beta-0n great-great granddaughter

If in your case you need more nuclei, feel free to add them as necessary. If you are not subject to beta-n decays, no need to delete those portions, just set the Pn and P2n probabilities to zero and they will not be filled.

The WriteToFile function will write all the histograms to a .root file for convenience.

To fit the data, use either DecayFit.C (which includes only a parent, daughter, and 2 contaminant components) or DecayFitwBetaN.C (which has components for the beta-n and beta-2n daughters plus granddaughters). Much of the fitting code was adopted from Ryan Dunlop (see github.com/r3dunlop).
