# Research questions

1) What if we run a BLAST for each gene to get a bunch of organisms, whatever they are?
- Ended up using the Ensembl homology API for each gene
2) Clean up the results from 2 by further sub-dividing the categories into sub-sets of organisms that are 
evolutionarily related based on phylogeny (Ensembl tree is a good resource for this and I am sure there are other for 
validation as well). This will probably increase the likelihood of us finding orthologs since, based on the phylogeny 
structure used, they should already be orthologs; By not limiting the study to mammals, reptiles, and cartilagenous 
fish, it might help the validity of the study. I think using only two organisms per category (based on my current 
understanding) would limit validity.
- This is related to 1. Ended up using the homology API, which returned numerous orthologs for each gene (the main 
species of interest was homo sapiens).
3) There are several genes that showcase many interspaced aligned regions. For example, Prg4 does not seem to showcase 
long regions of similar nucleotides. Rather, it shows many small regions that seem to align, with long break in between.
Why is this the case? What is the function of Prg4? Does it make sense for it to have high mutation rates such that we 
observe these inter-spaced regions, which can suggest inconsistency between organisms.
4)