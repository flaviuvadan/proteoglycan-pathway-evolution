import csv
import json
import os
import time

import requests

from src import exceptions


class Collector:
    """ Collector is responsible for collecting taxonomic information for all the unique organisms in data/organisms """

    GENE_ID_IDX = 1
    GENE_NAME_IDX = 0

    FREQ_KEY = "freq"
    GENE_IDS_KEY = "genes"

    def __init__(self):
        """
        Constructor
        """
        self.parent_path = self._get_parent_path()
        self.gene_file_path = self._get_gene_file_path()
        self.gene_ids = self._load_genes()
        # keep track of what genes an organism has and what organisms a gene appears in
        self.organisms, self.genes = self._load_orgs_and_genes()

    def _get_parent_path(self):
        """ Constructs the parent path of the folder where the files of this class are held """
        return os.path.join(os.getcwd(), "src", "data", "genes")

    def _get_gene_file_path(self):
        """ Constructs and returns the gene file paths """
        return os.path.join(os.getcwd(), "src", "data", "genes", "genes.txt")

    def _load_genes(self):
        """ Loads the genes into the class gene_ids """
        with open(self.gene_file_path, 'r') as gene_file:
            csv_reader = csv.reader(gene_file, delimiter=',')
            for gene in csv_reader:
                yield (gene[self.GENE_NAME_IDX], gene[self.GENE_ID_IDX])

    def _get_organisms_file_path(self, gene_name, gene_id):
        """ Builds the file path to the organisms file of the given gene name and ID """
        return os.path.join(os.getcwd(), "src", "data", "organisms", "{}_{}.txt".format(gene_name, gene_id))

    def _load_orgs_and_genes(self):
        """ Responsible for loading all the organisms curated for the list of genes of interest. Builds and returns two
        maps - one mapping each organism to its genes and one mapping each gene to its organisms """
        organisms = {}
        genes = {}
        for gene in self.gene_ids:
            org_file_path = self._get_organisms_file_path(gene[self.GENE_NAME_IDX], gene[self.GENE_ID_IDX])
            with open(org_file_path, "r") as orgs:
                org = orgs.read().splitlines()
                genes[gene[self.GENE_NAME_IDX]] = {}
                # we only care about unique organisms
                for o in org:
                    if not o.startswith(">"):
                        continue
                    clean_o = o.replace(">", "", 1).replace("_", " ").title()
                    # I hate to do this but there's a special case for Canis Familiaris
                    # EBI does not recognize it but it does recognize Canis Lupus (Canis Lupus Familiaris)
                    if "Canis Familiaris" in clean_o:
                        clean_o = "Canis lupus"
                    if not organisms.get(clean_o):
                        organisms[clean_o] = {self.FREQ_KEY: 1, self.GENE_IDS_KEY: [gene]}
                    else:
                        organisms[clean_o][self.FREQ_KEY] = organisms[clean_o][self.FREQ_KEY] + 1
                        organisms[clean_o][self.GENE_IDS_KEY].append(gene)
                    genes[gene[self.GENE_NAME_IDX]][clean_o] = 1
        return organisms, genes

    def get_gene_frequencies(self):
        """ Builds a file containing the gene frequencies of each organism """
        path = os.path.join(self.parent_path, "gene_frequencies.txt")
        with open(path, "w") as freqs:
            freqs.write("Organism,Gene Frequency\n")
            for org, data in self.organisms.items():
                freqs.write("{},{}\n".format(org, data.get(self.FREQ_KEY)))

    def get_genes_organisms(self):
        """ Builds a file containing which organisms a gene appears in """
        path = os.path.join(self.parent_path, "genes_organisms.txt")
        with open(path, "w") as f:
            f.write("Gene,Organisms\n")
            for gene in self.genes.keys():
                f.write("{},{}".format(gene, "/".join(self.genes.get(gene).keys()) + "\n"))

    def get_organisms_genes(self):
        """ Builds a file containing the genes exhibited by each organism """
        path = os.path.join(self.parent_path, "organisms_genes.txt")
        with open(path, "w") as freqs:
            freqs.write("Organism,Genes\n")
            for org, data in self.organisms.items():
                genes = ""
                for gene in data.get(self.GENE_IDS_KEY):
                    genes = genes + "{} ".format(gene[self.GENE_NAME_IDX])
                freqs.write("{},{}\n".format(org, genes))

    def collect(self):
        """ Collects taxa information from EBI (https://www.ebi.ac.uk/ena/browse/taxonomy-service) """
        path = os.path.join(os.getcwd(), "src", "data", "phylogeny", "taxa_info.txt")
        with open(path, "w") as f:
            f.write("Kingdom,Subkingdom,Phylum,Clade,Subphylum,Clade,Class,Subclass,Superorder,Order,Suborder," +
                    "Subsuborder,Family,Genus,Species\n")
            for org in self.organisms.keys():
                url = self._construct_ebi_taxon_url(org)
                req = requests.get(url)
                if req.status_code == 404:
                    print("Failed to find taxon information for organism: {}".format(org))
                elif req.status_code != 200:
                    message = "\nError calling EBI\n" + \
                              "expected status code 200, found: {}\n" + \
                              "address: {}\n" + \
                              "org: {}".format(req.status_code, url, org)
                    raise exceptions.SpeciesRequestException(message)
                parsed_req = json.loads(req.text)[0].get("lineage")  # single JSON object typically
                f.write(self._format_organism_info(org, parsed_req))
                # do not bombard EBI with 100s of API calls
                time.sleep(0.200)

    def _construct_ebi_taxon_url(self, org):
        """ Constructs and returns the EBI url """
        split = self._parse_organism_name(org)
        return "http://www.ebi.ac.uk/ena/data/taxonomy/v1/taxon/scientific-name/{}%20{}".format(split[0], split[1])

    def _parse_organism_name(self, org):
        """ Returns the split organism name as (genus, species) """
        split_org = org.split()
        return split_org[0], split_org[1]

    def _format_organism_info(self, org, info):
        """ Creates a nicely formatted string to add to organisms.txt """
        species_idx = 1
        split_org = self._parse_organism_name(org)
        return ",".join(info.replace(";", "").split()) + ",{}\n".format(split_org[species_idx])


if __name__ == "__main__":
    collector = Collector()
    collector.collect()
    collector.get_organisms_genes()
    collector.get_gene_frequencies()
    collector.get_genes_organisms()
