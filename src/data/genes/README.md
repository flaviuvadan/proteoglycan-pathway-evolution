# BLAST gene trees checks

Contains the gene tree checks against BLAST. 

##### Note

Some trees contains more than 20 organisms that do not have a specific gene. Since NCBI BLAST only allows a maximum of 20 organisms to be checked, if two organisms are part of a clade (speciated and show up as two children of a binary node), only one of them was used. This should be ok because the organisms are likely to be similar since they are part of the same clade and should only differ slightly from an evolutionary standpoint.
