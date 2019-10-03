import csv
import requests
import json
import os

from src.orthologs import exceptions


class Collector:
    """ Calls Ensembl APIs to collect orthologue information for a collection of genes """

    GENE_ID_IDX = 1

    def __init__(self, gene_filename):
        """
        Collector
        :param gene_filename: name of the file containing gene names
        """
        self.gene_filename = gene_filename
        self.gene_ids = self.load()

    def load(self):
        """ Loads the genes in the class gene_filename """
        with open(self.gene_filename, 'r') as gene_file:
            csv_reader = csv.reader(gene_file, delimiter=',')
            for gene in csv_reader:
                yield gene[self.GENE_ID_IDX]

    def construct_ensembl_url(self, gene_id):
        """
        Constructs the Ensembl URL to get homology information for a gene
        :param gene_id: Ensembl ID of the gene
        :return: Ensembl URL
        """
        return "https://rest.ensembl.org/homology/id/{}?type=orthologues;sequence=dna;cigar_line=0".format(gene_id)

    def collect(self):
        """ Launches the collection of the orthologues """
        headers = {"content-type": "application/json"}
        for gene_id in self.gene_ids:
            with open("{}.txt".format(gene_id), 'w') as out:
                url = self.construct_ensembl_url(gene_id)
                req = requests.get(url=url, headers=headers)
                if req.status_code != 200:
                    message = "Error calling Ensembl, expected status code 200, found: {}".format(req.status_code)
                    raise exceptions.OrthologRequestException(message)
                out.write(req.text)
                self.format(out.name)
                break

    def format(self, gene_filename):
        """ Invokes the json_formatter script with the given gene_filename """
        print("invoking json_formatter")
        os.system("bash json_formatter.sh {}".format(gene_filename))


if __name__ == "__main__":
    collector = Collector("genes.txt")
    collector.collect()