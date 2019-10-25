import csv
import json
import os
import time

import requests

from src import exceptions


class Collector:
    """ Collector is responsible for collecting taxonomic information for all the unique organisms in data/organisms"""

    GENE_ID_IDX = 1
    GENE_NAME_IDX = 0

    FREQ_KEY = "freq"
    GENE_IDS_KEY = "genes"

    def __init__(self, gene_file_path):
        """
        Constructor
        :param gene_file_path: name of the file containing gene names
        """
        self.gene_file_path = gene_file_path
        self.gene_ids = self.load_genes()
        self.organisms = self.load_organisms()

    def load_genes(self):
        """ Loads the genes into the class gene_ids """
        with open(self.gene_file_path, 'r') as gene_file:
            csv_reader = csv.reader(gene_file, delimiter=',')
            for gene in csv_reader:
                yield (gene[self.GENE_NAME_IDX], gene[self.GENE_ID_IDX])

    def get_organisms_file_path(self, gene_name, gene_id):
        """ Builds the file path to the organisms file of the given gene name and ID """
        return os.path.join(os.pardir, "data", "organisms", "{}_{}.txt".format(gene_name, gene_id))

    def load_organisms(self):
        """ Responsible for loading all the organisms curated for the list of genes of interest """
        organisms = {}
        for gene in self.gene_ids:
            org_file_path = self.get_organisms_file_path(gene[self.GENE_NAME_IDX], gene[self.GENE_ID_IDX])
            with open(org_file_path, "r") as orgs:
                org = orgs.read().splitlines()
                # we only care about unique organisms
                for o in org:
                    if not o.startswith(">"):
                        continue
                    clean_o = o.replace(">", "", 1).replace("_", " ").title()
                    if not organisms.get(clean_o):
                        organisms[clean_o] = {self.FREQ_KEY: 1, self.GENE_IDS_KEY: [gene]}
                    else:
                        organisms[clean_o][self.FREQ_KEY] = organisms[clean_o][self.FREQ_KEY] + 1
                        organisms[clean_o][self.GENE_IDS_KEY].append(gene)
        return organisms

    def get_gene_frequencies(self):
        """ Builds a file containing the gene frequencies of each organism """
        with open("gene_frequencies.txt", "w") as freqs:
            freqs.write("Organism,Gene Frequency\n")
            for org, data in self.organisms.items():
                freqs.write("{},{}\n".format(org, data.get(self.FREQ_KEY)))

    def get_organisms_genes(self):
        """ Builds a file containing the genes exhibited by each organism """
        with open("organisms_genes.txt", "w") as freqs:
            freqs.write("Organism,Genes\n")
            for org, data in self.organisms.items():
                genes = ""
                for gene in data.get(self.GENE_IDS_KEY):
                    genes = genes + "{} ".format(gene[self.GENE_NAME_IDX])
                freqs.write("{},{}\n".format(org, genes))

    def collect(self):
        """ Collects taxa information from GBIF (https://www.gbif.org/developer/species) """
        with open("organisms.txt", "w") as f:
            f.write(self.format_organism_info(["Organism", "Kingdom", "Phylum", "Order", "Family", "Genus", "Species"]))
            for org in self.organisms.keys():
                url, params = self.construct_gbif_url(org)
                req = requests.get(url, params=params)
                if req.status_code != 200:
                    message = "Error calling GBIF, expected status code 200, found: {}".format(req.status_code)
                    raise exceptions.SpeciesRequestException(message)
                f.write(self.parse_results(org, req))
                # do not bombard GBIF with 100s of API calls
                time.sleep(0.300)

    def construct_gbif_url(self, species):
        """ Constructs and returns the GBIF url """
        return "http://api.gbif.org/v1/species", {"name": species}

    def parse_results(self, organism, results):
        """ Parses the results from a single API call for taxa information of an organism and returns a nicely
        formatted string """
        results.text.replace("false", "False")  # this is a problem with the API
        result = json.loads(results.text).get("results")
        try:  # it can happen the results are empty
            org = result[0]
            info = [organism, org.get("kingdom"), org.get("phylum"), org.get("order"), org.get("family"),
                    org.get("genus"),
                    org.get("species")]
            return self.format_organism_info(info)
        except IndexError:
            info = [organism]  # appending this will at least tell us which ones we could not find and we can fix them
            return self.format_organism_info(info)

    def format_organism_info(self, info):
        """ Creates a nicely formatted string to add to organisms.txt """
        final = ""
        for i, val in enumerate(info):
            if i == len(info) - 1:
                final = final + "{}\n".format(val)
            else:
                final = final + "{},".format(val)
        return final


if __name__ == "__main__":
    gene_file_path = os.path.join(os.pardir, "data", "genes.txt")
    collector = Collector(gene_file_path)
    collector.collect()
    collector.get_organisms_genes()
    collector.get_gene_frequencies()
