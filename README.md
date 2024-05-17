# LB1_Project
### ASSIGNMENT:
Create a classifier for Kunitz-domain protein sequences using HMM profile built with structural information.

## CONSTRUCTION AND COMPARISON OF THREE HMM PROFILES FOR THE ANNOTATION OF PROTEINS BASED ON THE PRESENCE OF KUNITZ-TYPE PROTEASE INHIBITOR DOMAIN
### Abstract
Kunitz-type protease inhibitor domain is present in a lot of proteins and it gives them important functions, it’s in fact involved in protease inhibition and other biological functions. 
It’s a very conserved and studied domain because it can have several potential applications in different fields such as medicine, agriculture and biotechnology. Here, the aim is to build three classifiers (Hidden Markov Model (HMM) profiles) starting from structural information (retrieved with different queries), that are able to annotate proteins based on the prediction of the presence of a Kunitz domain. 
Using the resulting profiles to build a sequence Logo can be seen the conservation of important residues for the function and the structure. 
All the profiles report high performances and good classification capacity, but after the optimization and the tests on SwissProt sequences, using Matthew’s correlation coefficient as main performance metrics, the HMM profile starting from the query on PDB with the PFAM code obtained the best value, even if some problems in the annotation of some sequences were identified, so the performance should be re-evaluate after checking the annotation of all the sequences in the test set. 
