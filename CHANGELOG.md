# Documentation of project progress

## January 7th, 2020
- added more significant organisms based on feedback
- changed the logic of gaps vs. illegal characters in dnds/visualizer
- run visualizer for all the significant organisms again

## Dec 29th, 2019
- implemented the rest of the dnds visualizer
- fixed the dnds package (inconsistent naming of organisms between alignments and dnds was causing dnds scores to look 
odd)
- implemented a distribution visualization for the dN/dS values of all the significant organisms

## Dec 28th, 2019
- implemented the similarity mapper (uses Jaccard's index for similarity)
- implemented tree visualization for tree similarities 
- reimplemented the dnds package, turns out the dnds library does not work properly, will have to use the 
`ape` package in R (has 35 citations), turn the dnds package into a parsing and visualizing package, and
construct a pipe that calls the parsing script, calls R on the resulting files, then visualizes the dnds results
computed by the R `ape` package

## Dec 27th, 2019
- add more taxonomic details to the significant organisms phylogenetic tree
- performed a comparison against BLAST of organisms that are reported to not have a gene
- implemented a parser for the JSON files of the BLAST comparison

## Dec 26th, 2019
- chose and added a subset of the organisms that represent the classes of bones/cartilage for the study
- added the phylogenetic mapper for significant organisms

## Dec 3, 2019
- added first presentation

## Nov 28, 2019
- fixed an inaccurate piece of information regarding gene functions
- adjusted gene function chart to have a tight bounding box
- adjusted the gene frequencies chart - sorted organisms based on gene frequency 

## Nov 18, 2019
- start building the dnds package

## Nov 17, 2019
- adjusted MSA visualization to start x-axis at 0
- adjusted grouping of organisms based on taxonomic information

## Nov 4, 2019
- fixed the problem of false negative matches when building phylogeny trees for all the genes

## Nov 3, 2019
- added a phylogeny mapper for building the trees
- reorganized some data files
- wrote script for phylogenetic tree visualizations

## Nov 2, 2019
- implemented a taxon mapper for generating a plot of frequency of organisms per clade
- implemented a gene function visualizer

## Nov 1, 2019
- implemented the MSA parser
- moved some results files around (plots moved to src/data)

## Oct 31, 2019
- started building the MSAs parser
- noticed the alignments had a problem caused by the ClustalW format (limited number of characters for the organism 
name). Performed MSAs again, with Kalign on EBI, and saved in FASTA format 

## Oct 24, 2019
- make taxa collector build organism information files (gene frequency and genes per organism)

## Oct 13, 2019
- used EBI's Kalign to perform multiple sequence alignments for all genes
- added a test parser script for writing MSAs as single lines (will not use yet as they are challenging to interpret if 
written on single lines)
- built the taxa collector and collected organisms' information

## Oct 3, 2019
- moved exceptions out of orthologs package
- create curator package to parse homology information

## Oct 2, 2019
- curated gene Ensembl IDs list
- added orthologs package
- implemented Collector for calls to Ensembl

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

#### Potential programs to write
- parse the list of genes to get orthologs of each one, including species (using an Ensembl API);
- use the species information from 1 to hit some other API (taxonomy API from NCBI?) to identify vertebrates, chondrichthyes, etc and group them accordingly;
- further refine 2 by sub-grouping to create sub-sets of organisms that are close, based on evolution trees, taxa, or something else;
- get the dn/ds for all pairs of organisms in the sub-sets (organism 1 vs organism 2 dn/ds, 1 -> 3, 1 -> 4… etc)
