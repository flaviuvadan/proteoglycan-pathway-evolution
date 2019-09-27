# Documentation of project progress (changelog)

## Sep 21, 2019
- took the gene list form Brian and stored it in `genes.txt`
- used a pipe to clean up the list from the format:

```
gene1, gene2, ...
```
to

```
gene1 

gene2
```
via 

```
echo "gene1, gene2" | tr ',\s' '\n'
``` 
- looked into setting up a BLAST server and found [this resource](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- the objective is to find how a given set of genes (G) relate to cartilage production in multiple organisms (O)
- O will contain:

1. Organisms that make bone
1. Organisms that only produce cartilage
1. Produce neither cartilage nor bone
1. Used to have bone but now make cartilage (note, this is different than *only* producing cartilage)

- see `genes.txt` for G
- explored the [Ensembl Species Tree](https://uswest.ensembl.org/info/about/speciestree.html) to identify what organisms to use in the project (see `organisms.txt`)
- found the [REST documentation of Ensembl](https://rest.ensembl.org/documentation/info/homology_ensemblgene)
- since the objective at the start of the project is to find genomes of species that host the genes of interest, all the project has to do at the start is find the Ensembl IDs of the genes of interest and use `GET homology/id/:id` to get orthologs (the response includes organisms, which can then be grouped into bone/cartilage/neither/used to have bone)

#### Questions
- what if we run a BLAST for each gene to get a bunch of organisms, whatever they are
- take the results from 1 and  intersect them with organisms known to make bone, cartilage only, etc
- clean up the results from 2 by further sub-dividing the categories into sub-sets of organisms that are evolutionarily related based on phylogeny (Ensembl tree is a good resource for this and I am sure there are other for validation as well). This will probably increase the likelihood of us finding orthologs since, based on the phylogeny structure used, they should already be orthologs;
By not limiting the study to mammals, reptiles, and cartilagenous fish, it might help the validity of the study. I think using only two organisms per category (based on my current understanding) would limit validity

#### Potential programs to write
- parse the list of genes to get orthologs of each one, including species (using an Ensembl API);
- use the species information from 1 to hit some other API (taxonomy API from NCBI?) to identify vertebrates, chondrichthyes, etc and group them accordingly;
- further refine 2 by sub-grouping to create sub-sets of organisms that are close, based on evolution trees, taxa, or something else;
- get the dn/ds for all pairs of organisms in the sub-sets (organism 1 vs organism 2 dn/ds, 1 -> 3, 1 -> 4… etc)